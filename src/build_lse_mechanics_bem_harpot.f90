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


subroutine build_lse_mechanics_bem_harpot(kf,kr)

  ! Fortran 2003 intrinsic module
  use iso_fortran_env

  ! fbem modules
  use fbem_data_structures
  use fbem_string_handling
  use fbem_numerical
  use fbem_shape_functions
  use fbem_geometry
  use fbem_symmetry
  use fbem_bem_general
  use fbem_bem_stapot2d
  use fbem_bem_stapot3d
  use fbem_bem_harpot2d
  use fbem_bem_harpot3d
  use fbem_harpot_incident_field

  ! Module of problem variables
  use problem_variables

  ! No implicit variables
  implicit none

  ! I/O variables
  integer                           :: kf
  integer                           :: kr
  ! Local variables
  real(kind=real64)                 :: omega
  integer                           :: kb_int, sb_int
  logical                           :: sb_int_reversion
  integer                           :: sp_int
  integer                           :: ke_int, se_int
  integer                           :: se_int_n_nodes
  integer                           :: kn_col, sn_col
  integer                           :: kb, sb
  integer                           :: ke, se
  integer                           :: kn, knc, sn
  integer                           :: ss1, ss2
  integer                           :: kc
  real(kind=real64), allocatable    :: xi_i(:)
  ! Dataset at integration points
  logical                            :: node_freeterm_added(n_nodes)       ! Vector of flags to know if free-term has been calculated and added
  ! Region properties
  real(kind=real64)                  :: rho, d1J
  complex(kind=real64)               :: c, k
  type(fbem_bem_harpot3d_parameters) :: p3d
  type(fbem_bem_harpot2d_parameters) :: p2d
  ! Incident wave variables
  complex(kind=real64), allocatable :: p_inc(:), Un_inc(:)
  ! Kernels for SBIE integration
  complex(kind=real64), allocatable :: hp(:), gp(:)
  complex(kind=real64), allocatable :: hm(:), gm(:)
  real(kind=real64)                 :: c_plus, c_minus
  ! Kernels for HBIE integration
  complex(kind=real64), allocatable :: mp(:), lp(:)
  complex(kind=real64), allocatable :: mm(:), lm(:)
  ! Multiplier for Dual Burton & Miller formulation
  real(kind=real64)                 :: alpha
  ! Associated with free-term calculation
  real(kind=real64), allocatable    :: pphi_i(:)
  real(kind=real64), allocatable    :: sphi_i(:)
  integer                           :: n_c_elements
  real(kind=real64), allocatable    :: n_set_at_gn(:,:), n_set_at_gn_reversed(:,:)
  real(kind=real64), allocatable    :: t_set_at_gn(:,:), t_set_at_gn_reversed(:,:)
  ! Assembling control variable
  logical                           :: assemble
  ! Writing
  character(len=fbem_fmtstr)              :: fmtstr
  integer                                 :: output_fileunit
  character(len=fbem_filename_max_length) :: tmp_filename
  character(len=fbem_string_max_length)   :: tmp_string

  ! Message
  if (verbose_level.ge.2) then
    write(fmtstr,*) '(1x,a6,1x,i',fbem_nchar_int(region(kr)%id),',1x,a27)'
    call fbem_trimall(fmtstr)
    write(output_unit,fmtstr) 'Region', region(kr)%id, '(BE region, inviscid fluid)'
  end if

  ! Frequency
  omega=frequency(kf)

  ! ¿?¿?¿
  ! PONER ESTO EN UNA FUNCION, OJO! ESTA TAMBIEN EN CAMPO INCIDENTE
  ! ¿?¿?¿?






  ! Save the region properties to local variables
  rho=region(kr)%property_r(1)
  c=region(kr)%property_c(4)
  ! Wavenumber
  k=omega/c
  ! 1/J=rho*omega**2
  d1J=rho*omega**2
  ! Calculate the region parameters
  select case (problem%n)
    case (2)
      call fbem_bem_harpot2d_calculate_parameters(rho,c,omega,p2d)
    case (3)
      call fbem_bem_harpot3d_calculate_parameters(rho,c,omega,p3d)
  end select






  ! ============================== !
  ! EXPORT WAVE PROPAGATION SPEEDS !
  ! ============================== !

  if (export_wsp) then
    output_fileunit=fbem_get_valid_unit()
    tmp_filename=trim(output_filename)//'.wsp'
    call fbem_trim2b(tmp_filename)
    if ((kf.eq.1).and.(kr.eq.1)) then
      open(unit=output_fileunit,file=trim(tmp_filename),action='write',recl=fbem_file_record_length)
      write(output_fileunit,'(a)'  ) '# Program      : multifebe'
      write(output_fileunit,'(a)'  ) '# Version      : 1.0'
      write(output_fileunit,'(a)'  ) '# File_format  : wsp'
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
      write(output_fileunit,'(a)'  ) '# 1-2      Region id and region type (1: inviscid fluid, 2: viscoelastic, 3: poroelastic)'
      if (frequency_units.eq.'f') then
        write(output_fileunit,'(a)'  ) '# 3        Frequency (Hz)'
      else
        write(output_fileunit,'(a)'  ) '# 3        Frequency (rad/s)'
      end if
      select case (complex_notation)
      case (1)
        write(output_fileunit,'(a)'  )   '# 4-5      If column 2 == 1 (inviscid fluid region): Re(c), Im(c)'
        write(output_fileunit,'(a)'  )   '# 4-7      If column 2 == 2 (viscoelastic region): Re(cp), Im(cp), Re(cs), Im(cs)'
        write(output_fileunit,'(a)'  )   '# 4-9      If column 2 == 3 (poroelastic region): Re(cp1), Im(cp1), Re(cp2), Im(cp2), Re(cs), Im(cs)'
      case (2)
        write(output_fileunit,'(a)'  )   '# 4-5      If column 2 == 1 (inviscid fluid region): Abs(c), Arg(c)'
        write(output_fileunit,'(a)'  )   '# 4-7      If column 2 == 2 (viscoelastic region): Abs(cp), Arg(cp), Abs(cs), Arg(cs)'
        write(output_fileunit,'(a)'  )   '# 4-9      If column 2 == 3 (poroelastic region): Abs(cp1), Arg(cp1), Abs(cp2), Arg(cp2), Abs(cs), Arg(cs)'
      end select
    else
      open(unit=output_fileunit,file=trim(tmp_filename),access='append',recl=fbem_file_record_length)
    end if
    write(fmtstr,*) '(2',fmt_integer,')'
    call fbem_trimall(fmtstr)
    write (output_fileunit,fmtstr,advance='no') region(kr)%id, 1
    write(fmtstr,*) '(1',fmt_real,')'
    call fbem_trimall(fmtstr)
    if (frequency_units.eq.'f') then
      write (output_fileunit,fmtstr,advance='no') omega*c_1_2pi
    else
      write (output_fileunit,fmtstr,advance='no') omega
    end if
    write(fmtstr,*) '(2',fmt_real,')'
    call fbem_trimall(fmtstr)
    select case (complex_notation)
      case (1)
        write(output_fileunit,fmtstr) abs(c), fbem_zarg(c)
      case (2)
        write(output_fileunit,fmtstr) dreal(c), dimag(c)
    end select
    close(unit=output_fileunit)
  end if

  ! ==================================== !
  ! CALCULATE AND ASSEMBLE BEM INTEGRALS !
  ! ==================================== !

  ! Message
  if (verbose_level.ge.2) then
    write(fmtstr,*) '(a23,1x,i',fbem_nchar_int(region(kr)%id),')'
    call fbem_trimall(fmtstr)
    call fbem_timestamp_message(output_unit,2)
    write(output_unit,fmtstr) 'START assembling region',region(kr)%id
  end if

  !
  ! Loop through the BOUNDARIES of the REGION for INTEGRATION
  !
  do kb_int=1,region(kr)%n_boundaries
    sb_int=region(kr)%boundary(kb_int)
    sb_int_reversion=region(kr)%boundary_reversion(kb_int)
    sp_int=boundary(sb_int)%part
    !$omp parallel do schedule (dynamic) default (shared) private (se_int,se_int_n_nodes)
    do ke_int=1,part(sp_int)%n_elements
      se_int=part(sp_int)%element(ke_int)
      se_int_n_nodes=element(se_int)%n_nodes
      call build_lse_mechanics_bem_harpot_element(omega,kr,sb_int,sb_int_reversion,se_int,se_int_n_nodes,p2d,p3d)
    end do
    !$omp end parallel do
  end do

  ! Message
  if (verbose_level.ge.2) then
    write(fmtstr,*) '(a21,1x,i',fbem_nchar_int(region(kr)%id),')'
    call fbem_trimall(fmtstr)
    call fbem_timestamp_message(output_unit,2)
    write(output_unit,fmtstr) 'END assembling region',region(kr)%id
  end if

  ! ================================= !
  ! CALCULATE AND ASSEMBLE FREE-TERMS !
  ! ================================= !


  ! Falta modificar esto para half-space f.s.

  ! Message
  if (verbose_level.ge.2) then
    write(fmtstr,*) '(1x,a52)'
    call fbem_trimall(fmtstr)
    write(output_unit,fmtstr) 'Calculating and assembling analytical free-terms ...'
  end if

  ! Initialize the free-term control variable
  node_freeterm_added=.false.

  !
  ! Loop through the BOUNDARIES of the REGION
  !
  do kb_int=1,region(kr)%n_boundaries
    ! BOUNDARY
    sb_int=region(kr)%boundary(kb_int)
    sb_int_reversion=region(kr)%boundary_reversion(kb_int)
    sp_int=boundary(sb_int)%part

    !
    ! Loop through the ELEMENTS of the BOUNDARY
    !
    do ke_int=1,part(sp_int)%n_elements
      ! ELEMENT
      se_int=part(sp_int)%element(ke_int)
      se_int_n_nodes=element(se_int)%n_nodes
      ! Allocate element-wise variables
      allocate (hp(se_int_n_nodes),gp(se_int_n_nodes))
      allocate (hm(se_int_n_nodes),gm(se_int_n_nodes))
      allocate (mp(se_int_n_nodes),lp(se_int_n_nodes))
      allocate (mm(se_int_n_nodes),lm(se_int_n_nodes))
      allocate (p_inc(se_int_n_nodes),Un_inc(se_int_n_nodes))
      allocate (xi_i(element(se_int)%n_dimension),pphi_i(se_int_n_nodes),sphi_i(se_int_n_nodes))
      ! Copy the incident field over the element
      select case (boundary(sb_int)%coupling)
        case (fbem_boundary_coupling_be,fbem_boundary_coupling_be_fe)
          p_inc =element(se_int)%incident_c(1,:,1)
          Un_inc=element(se_int)%incident_c(2,:,1)
        case (fbem_boundary_coupling_be_be,fbem_boundary_coupling_be_fe_be)
          if (sb_int_reversion) then
            p_inc =element(se_int)%incident_c(1,:,2)
            Un_inc=element(se_int)%incident_c(2,:,2)
          else
            p_inc =element(se_int)%incident_c(1,:,1)
            Un_inc=element(se_int)%incident_c(2,:,1)
          end if
      end select

      !
      ! Loop through the NODES of the ELEMENT
      !
      do kn_col=1,se_int_n_nodes
        ! COLLOCATION NODE
        sn_col=element(se_int)%node(kn_col)
        ! Initialize assemble flag
        assemble=.false.

        ! True for dual formulations (Burton & Miller and DBEM) when the collocation point for SBIE and HBIE is the same.
        if (node(sn_col)%dual_is_common) then

          ! ========================================= !
          ! SBIE & HBIE AT THE SAME COLLOCATION POINT !
          ! ========================================= !

          ! Initialize
          assemble=.true.
          hp=0.
          gp=0.
          mp=0.
          lp=0.
          hm=0.
          gm=0.
          mm=0.
          lm=0.
          ! Calculate the shape functions vector at xi_i
          xi_i=element(se_int)%xi_i_hbie(:,kn_col)
          pphi_i=fbem_phi_hybrid(element(se_int)%type_f1,element(se_int)%delta_f,xi_i)
          sphi_i=fbem_phi_hybrid(element(se_int)%type_f2,element(se_int)%delta_f,xi_i)
          ! Add free-term to h+
          hp(:)=hp(:)+0.5d0*pphi_i
          lp(:)=lp(:)-0.5d0*sphi_i
          ! Add the free-term matrix of the inverted elements to h- if the integration boundary is a crack-like boundary.
          if (boundary(sb_int)%class.eq.fbem_boundary_class_cracklike) then
            hm(:)=hm(:)+0.5d0*pphi_i
            lm(:)=lm(:)+0.5d0*sphi_i
          end if
          ! If dual member is fbem_dual_burton_miller, add beta*HBIE to h and g respectively, these are used to assembling.
          if (node(sn_col)%dual.eq.fbem_dual_burton_miller) then
            alpha=node(sn_col)%alpha
            hp=hp+alpha*c_im/c*mp
            gp=gp+alpha*c_im/c*lp
            ! If the integration boundary is a crack-like boundary.
            if (boundary(sb_int)%class.eq.fbem_boundary_class_cracklike) then
              hm=hm+alpha*c_im/c*mm
              gm=gm+alpha*c_im/c*lm
            end if
          end if

        else

          ! ====== !
          !  SBIE  !
          ! ====== !

          ! If the col. has SBIE formulation and has not been calculated, then the SBIE terms are calculated and assembled.
          if ((node(sn_col)%sbie.eq.fbem_sbie).and.(.not.node_freeterm_added(sn_col))) then
            ! Initialize
            assemble=.true.
            node_freeterm_added(sn_col)=.true.
            hp=0.
            gp=0.
            hm=0.
            gm=0.
            xi_i=element(se_int)%xi_i_sbie(:,kn_col)
            !
            ! If the collocation point is in an edge or a vertex, the free-term has to be calculated.
            !
            if (fbem_check_xi_on_element_boundary(element(se_int)%type_g,xi_i)) then
              !
              ! Count the number of elements connected to the node for the integration region
              !
              ! Elements connected directly
              n_c_elements=node(sn_col)%n_elements
              ! Elements connected through common nodes
              do kn=1,node(sn_col)%n_nodes
                ! Selected common node
                sn=node(sn_col)%node(kn)
                ! If a "be" node
                if (part(node(sn)%part(1))%type.eq.fbem_part_be_boundary) then
                  ! If the boundary of the node is in the integration region, then add the number of elements
                  ! of the common node to the counter
                  sb=part(node(sn)%part(1))%entity
                  do kb=1,region(kr)%n_boundaries
                    if (region(kr)%boundary(kb).eq.sb) then
                      n_c_elements=n_c_elements+node(sn)%n_elements
                    end if
                  end do
                end if
              end do
              ! If the collocation node belongs to any symmetry plane
              select case (node(sn_col)%n_symplanes)
                ! If it belongs to 1 symmetry plane
                case (1)
                  n_c_elements=2*n_c_elements
                ! If it belongs to 2 symmetry planes
                case (2)
                  n_c_elements=4*n_c_elements
              end select
              ! Check
              if ((problem%n.eq.2).and.(n_c_elements.gt.2)) then
                call fbem_error_message(error_unit,0,'node',node(sn_col)%id,'is connected to more than 2 "be" elements.')
              end if
              ! Allocate the normals and tangents temporary matrices
              allocate (n_set_at_gn(problem%n,n_c_elements))
              allocate (t_set_at_gn(problem%n,n_c_elements))
              allocate (n_set_at_gn_reversed(problem%n,n_c_elements))
              allocate (t_set_at_gn_reversed(problem%n,n_c_elements))
              !
              ! Loop through the elements connected directly to the collocation node
              !
              ! Initialize the counter
              n_c_elements=0
              do ke=1,node(sn_col)%n_elements
                ! Selected element
                se=node(sn_col)%element(ke)
                ! The index of the node in the selected element
                kn=node(sn_col)%element_node_iid(ke)
                ! Increment the counter
                n_c_elements=n_c_elements+1
                ! Copy the normal and the tangent with the appropiate sign
                if (sb_int_reversion.eqv.(.false.)) then
                  do kc=1,problem%n
                    n_set_at_gn         (kc,n_c_elements)= element(se)%n_gn(kc,kn)
                    t_set_at_gn         (kc,n_c_elements)= element(se)%tbp_gn(kc,kn)
                    n_set_at_gn_reversed(kc,n_c_elements)=-element(se)%n_gn(kc,kn)
                    t_set_at_gn_reversed(kc,n_c_elements)= element(se)%tbm_gn(kc,kn)
                  end do
                else
                  do kc=1,problem%n
                    n_set_at_gn         (kc,n_c_elements)=-element(se)%n_gn(kc,kn)
                    t_set_at_gn         (kc,n_c_elements)= element(se)%tbm_gn(kc,kn)
                    n_set_at_gn_reversed(kc,n_c_elements)= element(se)%n_gn(kc,kn)
                    t_set_at_gn_reversed(kc,n_c_elements)= element(se)%tbp_gn(kc,kn)
                  end do
                end if
              end do
              !
              ! Loop through common nodes
              !
              do kn=1,node(sn_col)%n_nodes
                ! Selected common node
                sn=node(sn_col)%node(kn)
                ! If a "be" node
                if (part(node(sn)%part(1))%type.eq.fbem_part_be_boundary) then
                  ! Copy the boundary of the selected common node
                  sb=part(node(sn)%part(1))%entity
                  ! Loop through the boundaries of the integration region
                  do kb=1,region(kr)%n_boundaries
                    ! If the boundary of the selected common node is in the integration region, then copy all the
                    ! normals and tangents of the elements of the common node
                    if (region(kr)%boundary(kb).eq.sb) then
                      ! Loop through the elements of the common node
                      do ke=1,node(sn)%n_elements
                        ! Selected element
                        se=node(sn)%element(ke)
                        ! The index of the common node in the selected element
                        knc=node(sn)%element_node_iid(ke)
                        ! Increment the counter
                        n_c_elements=n_c_elements+1
                        ! Copy the normal and the tangent with the appropiate sign
                        if (region(kr)%boundary_reversion(kb).eqv.(.false.)) then
                          do kc=1,problem%n
                            n_set_at_gn         (kc,n_c_elements)= element(se)%n_gn(kc,knc)
                            t_set_at_gn         (kc,n_c_elements)= element(se)%tbp_gn(kc,knc)
                            n_set_at_gn_reversed(kc,n_c_elements)=-element(se)%n_gn(kc,knc)
                            t_set_at_gn_reversed(kc,n_c_elements)= element(se)%tbm_gn(kc,knc)
                          end do
                        else
                          do kc=1,problem%n
                            n_set_at_gn         (kc,n_c_elements)=-element(se)%n_gn(kc,knc)
                            t_set_at_gn         (kc,n_c_elements)= element(se)%tbm_gn(kc,knc)
                            n_set_at_gn_reversed(kc,n_c_elements)= element(se)%n_gn(kc,knc)
                            t_set_at_gn_reversed(kc,n_c_elements)= element(se)%tbp_gn(kc,knc)
                          end do
                        end if
                      end do
                    end if
                  end do
                end if
              end do
              !
              ! If the collocation node belongs to any symmetry plane, it is necessary to build the normals and
              ! tangents of symmetrical elements.
              !
              select case (node(sn_col)%n_symplanes)
                !
                ! If it belongs to 1 symmetry plane (it can happen in 2D and 3D)
                !
                case (1)
                  ! Symmetry plane of the collocation node
                  ss1=node(sn_col)%symplane(1)
                  ! Loop through the original elements
                  do ke=1,n_c_elements
                    ! Reflect the normal and the tangent with reversed orientation of each root element.
                    do kc=1,problem%n
                      n_set_at_gn         (kc,ke+n_c_elements)=symplane_m(kc,ss1)*n_set_at_gn         (kc,ke)
                      t_set_at_gn         (kc,ke+n_c_elements)=symplane_m(kc,ss1)*t_set_at_gn_reversed(kc,ke)
                      n_set_at_gn_reversed(kc,ke+n_c_elements)=symplane_m(kc,ss1)*n_set_at_gn_reversed(kc,ke)
                      t_set_at_gn_reversed(kc,ke+n_c_elements)=symplane_m(kc,ss1)*t_set_at_gn        (kc,ke)
                    end do
                  end do
                  ! Save the total number of elements
                  n_c_elements=2*n_c_elements
                !
                ! If it belongs to 2 symmetry planes (it can happen only in 3D)
                !
                case (2)
                  ! Symmetry planes of the collocation node
                  ss1=node(sn_col)%symplane(1)
                  ss2=node(sn_col)%symplane(2)
                  ! Loop through the original elements
                  do ke=1,n_c_elements
                    ! Reflect the normal and the tangent with reversed orientation of each root element with
                    ! respect to the first symmetry plane.
                    do kc=1,3
                      n_set_at_gn         (kc,ke+n_c_elements)=symplane_m(kc,ss1)*n_set_at_gn         (kc,ke)
                      t_set_at_gn         (kc,ke+n_c_elements)=symplane_m(kc,ss1)*t_set_at_gn_reversed(kc,ke)
                      n_set_at_gn_reversed(kc,ke+n_c_elements)=symplane_m(kc,ss1)*n_set_at_gn_reversed(kc,ke)
                      t_set_at_gn_reversed(kc,ke+n_c_elements)=symplane_m(kc,ss1)*t_set_at_gn         (kc,ke)
                    end do
                    ! Reflect the normal and the tangent of each root element with respect to the first and
                    ! the second symmetry planes.
                    do kc=1,3
                      n_set_at_gn         (kc,ke+2*n_c_elements)=symplane_m(kc,ss1)*symplane_m(kc,ss2)*n_set_at_gn         (kc,ke)
                      t_set_at_gn         (kc,ke+2*n_c_elements)=symplane_m(kc,ss1)*symplane_m(kc,ss2)*t_set_at_gn         (kc,ke)
                      n_set_at_gn_reversed(kc,ke+2*n_c_elements)=symplane_m(kc,ss1)*symplane_m(kc,ss2)*n_set_at_gn_reversed(kc,ke)
                      t_set_at_gn_reversed(kc,ke+2*n_c_elements)=symplane_m(kc,ss1)*symplane_m(kc,ss2)*t_set_at_gn_reversed(kc,ke)
                    end do
                    ! Reflect the normal and the tangent with reversed orientation of each root element with
                    ! respect to the second symmetry plane.
                    do kc=1,3
                      n_set_at_gn         (kc,ke+3*n_c_elements)=symplane_m(kc,ss2)*n_set_at_gn         (kc,ke)
                      t_set_at_gn         (kc,ke+3*n_c_elements)=symplane_m(kc,ss2)*t_set_at_gn_reversed(kc,ke)
                      n_set_at_gn_reversed(kc,ke+3*n_c_elements)=symplane_m(kc,ss2)*n_set_at_gn_reversed(kc,ke)
                      t_set_at_gn_reversed(kc,ke+3*n_c_elements)=symplane_m(kc,ss2)*t_set_at_gn         (kc,ke)
                    end do
                  end do
                  ! Save the total number of elements
                  n_c_elements=4*n_c_elements
              end select
              !
              ! Depending on the problem dimension, calculate the free-term of h+ and h-
              !
              select case (problem%n)
                case (2)
                  call fbem_bem_pot2d_sbie_freeterm(n_c_elements,n_set_at_gn,         t_set_at_gn,         geometric_tolerance,c_plus)
                  call fbem_bem_pot2d_sbie_freeterm(n_c_elements,n_set_at_gn_reversed,t_set_at_gn_reversed,geometric_tolerance,c_minus)
                case (3)
                  call fbem_bem_pot3d_sbie_freeterm(n_c_elements,n_set_at_gn,         t_set_at_gn,         geometric_tolerance,c_plus)
                  call fbem_bem_pot3d_sbie_freeterm(n_c_elements,n_set_at_gn_reversed,t_set_at_gn_reversed,geometric_tolerance,c_minus)
              end select
              ! Deallocate temporary variables
              deallocate (n_set_at_gn,t_set_at_gn)
              deallocate (n_set_at_gn_reversed,t_set_at_gn_reversed)
            !
            ! If the node is not in an edge or vertex of the element
            !
            else
              c_plus =0.5d0
              c_minus=0.5d0
            end if
            !
            ! Add the free-term matrix to h+ in the collocation node
            !
            hp(kn_col)=hp(kn_col)+c_plus
            ! Add the free-term matrix of the inverted elements to h- if the integration boundary is a crack-like boundary
            if (boundary(sb_int)%class.eq.fbem_boundary_class_cracklike) then
              hm(kn_col)=hm(kn_col)+c_minus
            end if
          end if

          ! ======== !
          ! SBIE MCA !
          ! ======== !

          ! If the collocation node has SBIE MCA formulation, then the SBIE MCA kernels of the integration element have to be
          ! calculated.
          if (node(sn_col)%sbie.eq.fbem_sbie_mca) then
            ! Initialize
            assemble=.true.
            hp=0.
            gp=0.
            hm=0.
            gm=0.
            ! Calculate the shape functions vector at xi_i
            xi_i=element(se_int)%xi_i_sbie_mca(:,kn_col)
            pphi_i=fbem_phi_hybrid(element(se_int)%type_f1,element(se_int)%delta_f,xi_i)
            ! Add free-term to h+
            hp(:)=hp(:)+0.5d0*pphi_i
            ! Add the free-term matrix of the inverted elements to h- if the integration boundary is a crack-like boundary.
            if (boundary(sb_int)%class.eq.fbem_boundary_class_cracklike) then
              hm(:)=hm(:)+0.5d0*pphi_i
            end if
          end if

          ! ====== !
          !  HBIE  !
          ! ====== !

          ! If the collocation node has HBIE (MCA) formulation, then the HBIE equation has to be integrated.
          if (node(sn_col)%hbie.eq.fbem_hbie) then
            ! Initialize
            assemble=.true.
            mp=0.
            lp=0.
            mm=0.
            lm=0.
            ! Calculate the shape functions vector at xi_i
            xi_i=element(se_int)%xi_i_hbie(:,kn_col)
            sphi_i=fbem_phi_hybrid(element(se_int)%type_f2,element(se_int)%delta_f,xi_i)
            ! Add free-term to h+
            lp(:)=lp(:)-0.5d0*sphi_i
            ! Add the free-term matrix of the inverted elements to h- if the integration boundary is a crack-like boundary.
            if (boundary(sb_int)%class.eq.fbem_boundary_class_cracklike) then
              lm(:)=lm(:)+0.5d0*sphi_i
            end if
            ! SBIE = HBIE or SBIE + beta*HBIE IF REQUIRED
            ! If dual member is 0, and the HBIE has been calculated, is because only the HBIE is used. Copy the
            ! m and l kernels to h and g respectively, because these are used to assembling.
            if (node(sn_col)%dual.eq.0) then
              hp=mp
              gp=lp
              ! If the integration boundary is a crack-like boundary.
              if (boundary(sb_int)%class.eq.fbem_boundary_class_cracklike) then
                hm=mm
                gm=lm
              end if
            end if
            ! If dual member is fbem_dual_burton_miller, add beta*HBIE to h and g respectively, these are used to assembling.
            if (node(sn_col)%dual.eq.fbem_dual_burton_miller) then
              alpha=node(sn_col)%alpha
              hp=hp+alpha*c_im/c*mp
              gp=gp+alpha*c_im/c*lp
              ! If the integration boundary is a crack-like boundary.
              if (boundary(sb_int)%class.eq.fbem_boundary_class_cracklike) then
                hm=hm+alpha*c_im/c*mm
                gm=gm+alpha*c_im/c*lm
              end if
            end if
          end if

        end if

        ! ======== !
        ! ASSEMBLE !
        ! ======== !

        ! Note: The flux variable is the normal displacement Un=1/(rho*omega**2)*dp/dn.
        ! The collocation establishes the equations (rows).
        ! The integration establishes the variables (columns).
        if (assemble) then
          select case (boundary(sb_int)%coupling)
            case (fbem_boundary_coupling_be,fbem_boundary_coupling_be_fe)
              select case (boundary(sb_int)%class)
                case (fbem_boundary_class_ordinary)
                  gp=gp*d1J
                  call assemble_bem_harpot_equation(omega,kr,sb_int,sb_int_reversion,se_int,se_int_n_nodes,sn_col,1,hp,gp,hm,gm,p_inc,Un_inc)
                case (fbem_boundary_class_cracklike)
                  gp=gp*d1J
                  gm=gm*d1J
                  lp=lp*d1J
                  lm=lm*d1J
                  call assemble_bem_harpot_equation(omega,kr,sb_int,sb_int_reversion,se_int,se_int_n_nodes,sn_col,1,hp,gp,hm,gm,p_inc,Un_inc)
                  call assemble_bem_harpot_equation(omega,kr,sb_int,sb_int_reversion,se_int,se_int_n_nodes,sn_col,2,mp,lp,mm,lm,p_inc,Un_inc)
              end select
            case (fbem_boundary_coupling_be_be,fbem_boundary_coupling_be_fe_be)
              gp=gp*d1J
              if (sb_int_reversion) then
                call assemble_bem_harpot_equation(omega,kr,sb_int,sb_int_reversion,se_int,se_int_n_nodes,sn_col,2,hp,gp,hm,gm,p_inc,Un_inc)
              else
                call assemble_bem_harpot_equation(omega,kr,sb_int,sb_int_reversion,se_int,se_int_n_nodes,sn_col,1,hp,gp,hm,gm,p_inc,Un_inc)
              end if
          end select
        end if

      end do ! Loop through the NODES of the ELEMENT

      ! Deallocate element-wise data structures
      deallocate (hp,gp,hm,gm)
      deallocate (mp,lp,mm,lm)
      deallocate (p_inc,Un_inc)
      deallocate (xi_i,pphi_i,sphi_i)

    end do ! Loop through the ELEMENTS of the BOUNDARY

  end do ! Loop through the BOUNDARIES of the REGION

  ! Ending message
  if (verbose_level.ge.2) write(output_unit,'(1x,a)') 'done.'

end subroutine build_lse_mechanics_bem_harpot

subroutine build_lse_mechanics_bem_harpot_element(omega,kr,sb_int,sb_int_reversion,se_int,se_int_n_nodes,p2d,p3d)

  ! Fortran 2003 intrinsic module
  use iso_fortran_env

  ! fbem modules
  use fbem_data_structures
  use fbem_string_handling
  use fbem_numerical
  use fbem_symmetry
  use fbem_bem_general
  use fbem_bem_harpot2d
  use fbem_bem_harpot3d

  ! Module of problem variables
  use problem_variables

  ! No implicit variables
  implicit none

  ! I/O variables
  real(kind=real64)                  :: omega
  integer                            :: kr
  integer                            :: sb_int
  logical                            :: sb_int_reversion
  integer                            :: se_int
  integer                            :: se_int_n_nodes
  type(fbem_bem_harpot3d_parameters) :: p3d
  type(fbem_bem_harpot2d_parameters) :: p2d
  ! Local variables
  integer                :: ks
  type(fbem_bem_element) :: se_int_data
  logical                :: se_int_reversion
  integer                :: kb_col, sb_col
  logical                :: sb_col_reversion
  integer                :: sp_col
  integer                :: ke_col, se_col
  integer                :: se_col_n_nodes
  integer                :: kn_col, sn_col
  integer                :: kn
  real(kind=real64)      :: x_i(problem%n), n_i(problem%n)
  ! Dataset at integration points
  logical                :: node_collocated(n_nodes)
  ! Region properties
  real(kind=real64)      :: rho
  complex(kind=real64)   :: c, k
  real(kind=real64)      :: d1J
  ! Incident wave variables
  complex(kind=real64)   :: p_inc(se_int_n_nodes), Un_inc(se_int_n_nodes)
  ! Kernels for SBIE integration
  complex(kind=real64)   :: h (se_int_n_nodes), g (se_int_n_nodes)
  complex(kind=real64)   :: hp(se_int_n_nodes), gp(se_int_n_nodes)
  complex(kind=real64)   :: hm(se_int_n_nodes), gm(se_int_n_nodes)
  ! Kernels for HBIE integration
  complex(kind=real64)   :: m (se_int_n_nodes), l (se_int_n_nodes)
  complex(kind=real64)   :: mp(se_int_n_nodes), lp(se_int_n_nodes)
  complex(kind=real64)   :: mm(se_int_n_nodes), lm(se_int_n_nodes)
  ! Multiplier for Dual Burton & Miller formulation
  real(kind=real64)      :: alpha
  ! Associated with symmetry
  real(kind=real64)      :: symconf_m(problem%n), symconf_t(problem%n), symconf_r(problem%n), symconf_s
  logical                :: reversed
  ! Assembling control variable
  logical                :: assemble

  ! Save the region properties to local variables
  rho=region(kr)%property_r(1)
  select case (problem%n)
    case (2)
      c=p2d%c
    case (3)
      c=p3d%c
  end select
  k=omega/c
  d1J=rho*omega**2

  ! Initialize calculation element
  call se_int_data%init
  se_int_data%gtype=element(se_int)%type_g
  se_int_data%d=element(se_int)%n_dimension
  se_int_data%n_gnodes=element(se_int)%n_nodes
  se_int_data%n=problem%n
  allocate (se_int_data%x(problem%n,se_int_n_nodes))
  se_int_data%x=element(se_int)%x_gn
  se_int_data%ptype=element(se_int)%type_f1
  se_int_data%ptype_delta=element(se_int)%delta_f
  se_int_data%n_pnodes=element(se_int)%n_nodes
  se_int_data%stype=element(se_int)%type_f2
  se_int_data%stype_delta=element(se_int)%delta_f
  se_int_data%n_snodes=element(se_int)%n_nodes
  se_int_data%cl=element(se_int)%csize
  se_int_data%gln_far=element(se_int)%n_phi
  allocate (se_int_data%bball_centre(problem%n))
  se_int_data%bball_centre=element(se_int)%bball_centre
  se_int_data%bball_radius=element(se_int)%bball_radius
  ! Copy the incident field over the element
  select case (boundary(sb_int)%coupling)
    case (fbem_boundary_coupling_be,fbem_boundary_coupling_be_fe)
      p_inc =element(se_int)%incident_c(1,:,1)
      Un_inc=element(se_int)%incident_c(2,:,1)
    case (fbem_boundary_coupling_be_be,fbem_boundary_coupling_be_fe_be)
      if (sb_int_reversion) then
        p_inc =element(se_int)%incident_c(1,:,2)
        Un_inc=element(se_int)%incident_c(2,:,2)
      else
        p_inc =element(se_int)%incident_c(1,:,1)
        Un_inc=element(se_int)%incident_c(2,:,1)
      end if
  end select

  !
  ! Loop through symmetrical elements
  !
  do ks=1,n_symelements
    ! SYMMETRY SETUP
    call fbem_symmetry_multipliers(ks,problem%n,n_symplanes,symplane_m,symplane_s,symplane_t,symplane_r,&
                                   symconf_m,symconf_s,symconf_t,symconf_r,reversed)
    ! Change of element orientation and coordinates due to symmetry
    se_int_reversion=sb_int_reversion.neqv.reversed
    do kn=1,se_int_n_nodes
      se_int_data%x(:,kn)=symconf_m*element(se_int)%x_gn(:,kn)
    end do
    ! Initialize precalculated datasets
    call se_int_data%init_precalculated_datasets(n_precalsets,precalset_gln)

    !
    ! Loop through the BOUNDARIES of the REGION for COLLOCATION
    !
    ! Initialize the collocation control variable used when a node has the SBIE formulation
    node_collocated=.false.
    do kb_col=1,region(kr)%n_boundaries
      sb_col=region(kr)%boundary(kb_col)
      sb_col_reversion=region(kr)%boundary_reversion(kb_col)
      sp_col=boundary(sb_col)%part

      !
      ! Loop through the ELEMENTS of the BOUNDARY for COLLOCATION
      !
      do ke_col=1,part(sp_col)%n_elements
        ! COLLOCATION ELEMENT
        se_col=part(sp_col)%element(ke_col)
        se_col_n_nodes=element(se_col)%n_nodes

        !
        ! Loop through the NODES of the ELEMENT for COLLOCATION
        !
        do kn_col=1,se_col_n_nodes
          ! COLLOCATION NODE
          sn_col=element(se_col)%node(kn_col)
          ! Initialize assemble flag
          assemble=.false.

          ! True for dual formulations (Burton & Miller and DBEM) when the collocation point for SBIE and HBIE is the same.
          if (node(sn_col)%dual_is_common) then

            ! ==========================================
            !  SBIE & HBIE AT THE SAME COLLOCATION POINT
            ! ==========================================

            ! Initialize
            assemble=.true.
            ! CALCULATE KERNELS
            x_i=element(se_col)%x_i_hbie(:,kn_col)
            n_i=element(se_col)%n_i_hbie(:,kn_col)
            if (sb_col_reversion) n_i=-n_i
            select case (problem%n)
              case (2)
                call fbem_bem_harpot2d_shbie_auto(se_int_data,se_int_reversion,x_i,n_i,p2d,qsi_parameters,qsi_ns_max,h,g,m,l)
              case (3)
                call fbem_bem_harpot3d_shbie_auto(se_int_data,se_int_reversion,x_i,n_i,p3d,qsi_parameters,qsi_ns_max,h,g,m,l)
            end select
            ! Additional kernels for half-space fundamental solution
            if (region(kr)%space.eq.fbem_half_space) then
              x_i(abs(region(kr)%halfspace_n))=2.d0*region(kr)%halfspace_x-x_i(abs(region(kr)%halfspace_n))
              n_i(abs(region(kr)%halfspace_n))=-n_i(abs(region(kr)%halfspace_n))
              select case (problem%n)
                case (2)
                  call fbem_bem_harpot2d_shbie_auto(se_int_data,se_int_reversion,x_i,n_i,p2d,qsi_parameters,qsi_ns_max,hp,gp,mp,lp)
                case (3)
                  call fbem_bem_harpot3d_shbie_auto(se_int_data,se_int_reversion,x_i,n_i,p3d,qsi_parameters,qsi_ns_max,hp,gp,mp,lp)
              end select
              select case (region(kr)%halfspace_bc)
                ! p=0
                case (0)
                  h=h-hp
                  g=g-gp
                  m=m-mp
                  l=l-lp
                ! Un=0
                case (1)
                  h=h+hp
                  g=g+gp
                  m=m+mp
                  l=l+lp
              end select
            end if
            ! BUILD KERNELS ACCORDING TO SYMMETRY
            if (ks.gt.1) then
              h(:)=symconf_s*h(:)
              g(:)=symconf_s*g(:)
              m(:)=symconf_s*m(:)
              l(:)=symconf_s*l(:)
            end if
            ! BUILD KERNELS WITH N+ AND N-
            hp=h
            gp=g
            mp=m
            lp=l
            ! If the integration boundary is a crack-like boundary, build N- kernels
            if (boundary(sb_int)%class.eq.fbem_boundary_class_cracklike) then
              hm=-h
              gm= g
              mm=-m
              lm= l
            end iF
            ! If dual member is fbem_dual_burton_miller, add beta*HBIE to h and g respectively, these are used to assembling.
            if (node(sn_col)%dual.eq.fbem_dual_burton_miller) then
              alpha=node(sn_col)%alpha
              hp=hp+alpha*c_im/c*mp
              gp=gp+alpha*c_im/c*lp
              ! If the integration boundary is a crack-like boundary.
              if (boundary(sb_int)%class.eq.fbem_boundary_class_cracklike) then
                hm=hm+alpha*c_im/c*mm
                gm=gm+alpha*c_im/c*lm
              end if
            end if

          else

            ! ======
            !  SBIE
            ! ======

            ! If the collocation node has SBIE formulation, then the SBIE kernels of the integration element have to be
            ! calculated.

            if (node(sn_col)%sbie.eq.fbem_sbie.and.(.not.node_collocated(sn_col))) then
              ! Initialize
              assemble=.true.
              node_collocated(sn_col)=.true.
              ! CALCULATE KERNELS
              x_i=element(se_col)%x_i_sbie(:,kn_col)
              select case (problem%n)
                case (2)
                  call fbem_bem_harpot2d_sbie_auto(se_int_data,se_int_reversion,x_i,p2d,qsi_parameters,qsi_ns_max,h,g)
                case (3)
                  call fbem_bem_harpot3d_sbie_auto(se_int_data,se_int_reversion,x_i,p3d,qsi_parameters,qsi_ns_max,h,g)
              end select
              ! Additional kernels for half-space fundamental solution
              if (region(kr)%space.eq.fbem_half_space) then
                x_i(abs(region(kr)%halfspace_n))=2.d0*region(kr)%halfspace_x-x_i(abs(region(kr)%halfspace_n))
                select case (problem%n)
                  case (2)
                    call fbem_bem_harpot2d_sbie_auto(se_int_data,se_int_reversion,x_i,p2d,qsi_parameters,qsi_ns_max,hp,gp)
                  case (3)
                    call fbem_bem_harpot3d_sbie_auto(se_int_data,se_int_reversion,x_i,p3d,qsi_parameters,qsi_ns_max,hp,gp)
                end select
                select case (region(kr)%halfspace_bc)
                  ! p=0
                  case (0)
                    h=h-hp
                    g=g-gp
                  ! Un=0
                  case (1)
                    h=h+hp
                    g=g+gp
                end select
              end if
              ! BUILD KERNELS ACCORDING TO SYMMETRY
              if (ks.gt.1) then
                h(:)=symconf_s*h(:)
                g(:)=symconf_s*g(:)
              end if
              ! BUILD KERNELS WITH N+ AND N-
              hp=h
              gp=g
              ! If the integration boundary is a crack-like boundary, build N- kernels
              if (boundary(sb_int)%class.eq.fbem_boundary_class_cracklike) then
                hm=-h
                gm= g
              end if
            end if

            ! ==========
            !  SBIE MCA
            ! ==========

            ! If the collocation node has SBIE MCA formulation, then the SBIE MCA kernels of the integration element have to be
            ! calculated.
            if (node(sn_col)%sbie.eq.fbem_sbie_mca) then
              ! Initialize
              assemble=.true.
              ! CALCULATE KERNELS
              x_i=element(se_col)%x_i_sbie_mca(:,kn_col)
              select case (problem%n)
                case (2)
                  call fbem_bem_harpot2d_sbie_auto(se_int_data,se_int_reversion,x_i,p2d,qsi_parameters,qsi_ns_max,h,g)
                case (3)
                  call fbem_bem_harpot3d_sbie_auto(se_int_data,se_int_reversion,x_i,p3d,qsi_parameters,qsi_ns_max,h,g)
              end select
              ! Additional kernels for half-space fundamental solution
              if (region(kr)%space.eq.fbem_half_space) then
                x_i(abs(region(kr)%halfspace_n))=2.d0*region(kr)%halfspace_x-x_i(abs(region(kr)%halfspace_n))
                select case (problem%n)
                  case (2)
                    call fbem_bem_harpot2d_sbie_auto(se_int_data,se_int_reversion,x_i,p2d,qsi_parameters,qsi_ns_max,hp,gp)
                  case (3)
                    call fbem_bem_harpot3d_sbie_auto(se_int_data,se_int_reversion,x_i,p3d,qsi_parameters,qsi_ns_max,hp,gp)
                end select
                select case (region(kr)%halfspace_bc)
                  ! p=0
                  case (0)
                    h=h-hp
                    g=g-gp
                  ! Un=0
                  case (1)
                    h=h+hp
                    g=g+gp
                end select
              end if
              ! BUILD KERNELS ACCORDING TO SYMMETRY
              if (ks.gt.1) then
                h(:)=symconf_s*h(:)
                g(:)=symconf_s*g(:)
              end if
              ! BUILD KERNELS WITH N+ AND N-
              hp=h
              gp=g
              ! If the integration boundary is a crack-like boundary, build N- kernels
              if (boundary(sb_int)%class.eq.fbem_boundary_class_cracklike) then
                hm=-h
                gm= g
              end if
            end if

            ! ======
            !  HBIE
            ! ======

            ! If the collocation node has HBIE (MCA) formulation, then the HBIE equation has to be integrated.
            if (node(sn_col)%hbie.eq.fbem_hbie) then
              ! Initialize
              assemble=.true.
              ! CALCULATE KERNELS
              x_i=element(se_col)%x_i_hbie(:,kn_col)
              n_i=element(se_col)%n_i_hbie(:,kn_col)
              if (sb_col_reversion) n_i=-n_i
              select case (problem%n)
                case (2)
                  call fbem_bem_harpot2d_hbie_auto(se_int_data,se_int_reversion,x_i,n_i,p2d,qsi_parameters,qsi_ns_max,m,l)
                case (3)
                  call fbem_bem_harpot3d_hbie_auto(se_int_data,se_int_reversion,x_i,n_i,p3d,qsi_parameters,qsi_ns_max,m,l)
              end select
              ! Additional kernels for half-space fundamental solution
              if (region(kr)%space.eq.fbem_half_space) then
                x_i(abs(region(kr)%halfspace_n))=2.d0*region(kr)%halfspace_x-x_i(abs(region(kr)%halfspace_n))
                n_i(abs(region(kr)%halfspace_n))=-n_i(abs(region(kr)%halfspace_n))
                select case (problem%n)
                  case (2)
                    call fbem_bem_harpot2d_hbie_auto(se_int_data,se_int_reversion,x_i,n_i,p2d,qsi_parameters,qsi_ns_max,mp,lp)
                  case (3)
                    call fbem_bem_harpot3d_hbie_auto(se_int_data,se_int_reversion,x_i,n_i,p3d,qsi_parameters,qsi_ns_max,mp,lp)
                end select
                select case (region(kr)%halfspace_bc)
                  ! p=0
                  case (0)
                    m=m-mp
                    l=l-lp
                  ! Un=0
                  case (1)
                    m=m+mp
                    l=l+lp
                end select
              end if
              ! BUILD KERNELS ACCORDING TO SYMMETRY
              if (ks.gt.1) then
                  m(:)=symconf_s*m(:)
                  l(:)=symconf_s*l(:)
              end if
              ! BUILD KERNELS WITH N+ AND N-
              mp=m
              lp=l
              ! If the integration boundary is a crack-like boundary, build N- kernels
              if (boundary(sb_int)%class.eq.fbem_boundary_class_cracklike) then
                mm=-m
                lm= l
              end if
              ! SBIE = HBIE or SBIE + beta*HBIE IF REQUIRED
              ! If dual member is 0, and the HBIE has been calculated, is because only the HBIE is used. Copy the
              ! m and l kernels to h and g respectively, because these are used to assembling.
              if (node(sn_col)%dual.eq.0) then
                hp=mp
                gp=lp
                ! If the integration boundary is a crack-like boundary.
                if (boundary(sb_int)%class.eq.fbem_boundary_class_cracklike) then
                  hm=mm
                  gm=lm
                end if
              end if
              ! If dual member is fbem_dual_burton_miller, add beta*HBIE to h and g respectively, these are used to assembling.
              if (node(sn_col)%dual.eq.fbem_dual_burton_miller) then
                alpha=node(sn_col)%alpha
                hp=hp+alpha*c_im/c*mp
                gp=gp+alpha*c_im/c*lp
                ! If the integration boundary is a crack-like boundary.
                if (boundary(sb_int)%class.eq.fbem_boundary_class_cracklike) then
                  hm=hm+alpha*c_im/c*mm
                  gm=gm+alpha*c_im/c*lm
                end if
              end if
            end if

          end if

          ! ======== !
          ! ASSEMBLE !
          ! ======== !

          ! Note: The flux variable is the normal displacement Un=1/(rho*omega**2)*dp/dn.
          ! The collocation establishes the equations (rows).
          ! The integration establishes the variables (columns).
          if (assemble) then
            !$omp critical
            select case (boundary(sb_col)%coupling)
              case (fbem_boundary_coupling_be,fbem_boundary_coupling_be_fe)
                select case (boundary(sb_col)%class)
                  case (fbem_boundary_class_ordinary)
                    gp=gp*d1J
                    call assemble_bem_harpot_equation(omega,kr,sb_int,sb_int_reversion,se_int,se_int_n_nodes,sn_col,1,hp,gp,hm,gm,p_inc,Un_inc)
                  case (fbem_boundary_class_cracklike)
                    gp=gp*d1J
                    gm=gm*d1J
                    lp=lp*d1J
                    lm=lm*d1J
                    call assemble_bem_harpot_equation(omega,kr,sb_int,sb_int_reversion,se_int,se_int_n_nodes,sn_col,1,hp,gp,hm,gm,p_inc,Un_inc)
                    call assemble_bem_harpot_equation(omega,kr,sb_int,sb_int_reversion,se_int,se_int_n_nodes,sn_col,2,mp,lp,mm,lm,p_inc,Un_inc)
                end select
              case (fbem_boundary_coupling_be_be,fbem_boundary_coupling_be_fe_be)
                gp=gp*d1J
                if (sb_col_reversion) then
                  call assemble_bem_harpot_equation(omega,kr,sb_int,sb_int_reversion,se_int,se_int_n_nodes,sn_col,2,hp,gp,hm,gm,p_inc,Un_inc)
                else
                  call assemble_bem_harpot_equation(omega,kr,sb_int,sb_int_reversion,se_int,se_int_n_nodes,sn_col,1,hp,gp,hm,gm,p_inc,Un_inc)
                end if
            end select
            !$omp end critical
          end if

        end do ! Loop through the NODES of the ELEMENT for COLLOCATION

      end do ! Loop through the ELEMENTS of the BOUNDARY for COLLOCATION

    end do ! Loop through the BOUNDARIES of the REGION for COLLOCATION

  end do ! Loop through SYMMETRICAL ELEMENTS

end subroutine build_lse_mechanics_bem_harpot_element
