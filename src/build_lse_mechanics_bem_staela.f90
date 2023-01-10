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

subroutine build_lse_mechanics_bem_staela(kr)

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
  use fbem_bem_staela2d
  use fbem_bem_staela3d

  ! Module of problem variables
  use problem_variables

  ! No implicit variables
  implicit none

  ! I/O variables
  integer                           :: kr
  ! Local variables
  integer                           :: il, ik
  integer                           :: kb_int, sb_int
  logical                           :: sb_int_reversion
  integer                           :: sp_int
  integer                           :: ke_int, se_int
  integer                           :: se_int_n_nodes, se_int_n_snodes
  integer                           :: kb, sb
  integer                           :: ke, se
  integer                           :: kn_col, sn_col, sn_fe
  integer                           :: kn_int, sn_int
  integer                           :: kn, knc, sn
  integer                           :: ss1, ss2
  integer                           :: kc
  real(kind=real64)                 :: x_i(problem%n)
  real(kind=real64), allocatable    :: xi_i(:)
  ! Dataset at integration points
  logical                           :: node_freeterm_added(n_nodes)       ! Vector of flags to know if free-term has been calculated and added
  ! Region properties
  real(kind=real64)                 :: mu, nu

  ! Kernels for SBIE integration
  real(kind=real64), allocatable    :: h (:,:,:), g (:,:,:)
  real(kind=real64), allocatable    :: hp(:,:,:), gp(:,:,:)
  real(kind=real64), allocatable    :: hm(:,:,:), gm(:,:,:)
  real(kind=real64)                 :: c_plus(problem%n,problem%n),c_minus(problem%n,problem%n)
  ! Kernels for HBIE integration
  real(kind=real64), allocatable    :: m (:,:,:), l (:,:,:)
  real(kind=real64), allocatable    :: mp(:,:,:), lp(:,:,:)
  real(kind=real64), allocatable    :: mm(:,:,:), lm(:,:,:)
  ! Multiplier for Dual Burton & Miller formulation
  real(kind=real64)                 :: alpha
  ! Associated with free-term calculation
  real(kind=real64), allocatable    :: pphi_i(:)
  real(kind=real64), allocatable    :: sphi_i(:)
  integer                           :: n_c_elements
  real(kind=real64), allocatable    :: n_set_at_gn(:,:), n_set_at_gn_reversed(:,:)
  real(kind=real64), allocatable    :: t_set_at_gn(:,:), t_set_at_gn_reversed(:,:)
  ! Assembling control variable
  integer                           :: row, col
  logical                           :: assemble
  ! Writing
  character(len=fbem_fmtstr)              :: fmtstr
  integer                                 :: output_fileunit
  character(len=fbem_filename_max_length) :: tmp_filename

  if (verbose_level.ge.2) then
    write(fmtstr,*) '(a26,1x,i',fbem_nchar_int(region(kr)%id),')'
    call fbem_trimall(fmtstr)
    call fbem_timestamp_message(output_unit,2)
    write(output_unit,fmtstr) 'START assembling BE region',region(kr)%id
  end if

  ! Copy material properties
  mu=region(kr)%property_r(2)
  nu=region(kr)%property_r(3)
  if ((problem%n.eq.2).and.(problem%subtype.eq.fbem_mechanics_plane_stress)) nu=nu/(1.d0+nu)

  ! ================================================================================================================================
  ! CALCULATE AND ASSEMBLE BEM INTEGRALS
  ! ================================================================================================================================

  ! -------------------------------- !
  ! INTEGRATION of REGION BOUNDARIES !
  ! -------------------------------- !

  do kb_int=1,region(kr)%n_boundaries
    sb_int=region(kr)%boundary(kb_int)
    sb_int_reversion=region(kr)%boundary_reversion(kb_int)
    sp_int=boundary(sb_int)%part
    !$omp parallel do schedule (dynamic) default (shared) private (se_int,se_int_n_nodes,fmtstr)
    do ke_int=1,part(sp_int)%n_elements
      se_int=part(sp_int)%element(ke_int)
      se_int_n_nodes=element(se_int)%n_nodes
      call build_lse_mechanics_bem_staela_element(kr,sb_int,sb_int_reversion,se_int,se_int_n_nodes,mu,nu)
    end do
    !$omp end parallel do
  end do

  ! ----------------------------------- !
  ! INTEGRATION of REGION BE BODY LOADS !
  ! ----------------------------------- !

  do kb_int=1,region(kr)%n_be_bodyloads
    sb_int=region(kr)%be_bodyload(kb_int)
    sp_int=be_bodyload(sb_int)%part
    !$omp parallel do schedule (dynamic) default (shared) private (se_int,se_int_n_nodes,fmtstr)
    do ke_int=1,part(sp_int)%n_elements
      se_int=part(sp_int)%element(ke_int)
      se_int_n_nodes=element(se_int)%n_nodes
      call build_lse_mechanics_bem_staela_bl(kr,sb_int,se_int,se_int_n_nodes,mu,nu)
    end do
    !$omp end parallel do
  end do

  ! ================================================================================================================================
  ! CALCULATE AND ASSEMBLE FREE-TERMS
  ! ================================================================================================================================

  ! Message
  if (verbose_level.ge.2) then
    write(fmtstr,*) '(1x,a41)'
    call fbem_trimall(fmtstr)
    write(output_unit,fmtstr) 'Calculating and assembling free-terms ...'
  end if

  ! Initialize the free-term control variable
  node_freeterm_added=.false.

  ! -------------------------------- !
  ! INTEGRATION of REGION BOUNDARIES !
  ! -------------------------------- !

  do kb_int=1,region(kr)%n_boundaries
    sb_int=region(kr)%boundary(kb_int)
    sb_int_reversion=region(kr)%boundary_reversion(kb_int)
    sp_int=boundary(sb_int)%part
    do ke_int=1,part(sp_int)%n_elements
      se_int=part(sp_int)%element(ke_int)
      se_int_n_nodes=element(se_int)%n_nodes

      allocate (hp(se_int_n_nodes,problem%n,problem%n),gp(se_int_n_nodes,problem%n,problem%n))
      allocate (hm(se_int_n_nodes,problem%n,problem%n),gm(se_int_n_nodes,problem%n,problem%n))
      allocate (mp(se_int_n_nodes,problem%n,problem%n),lp(se_int_n_nodes,problem%n,problem%n))
      allocate (mm(se_int_n_nodes,problem%n,problem%n),lm(se_int_n_nodes,problem%n,problem%n))
      allocate (xi_i(element(se_int)%n_dimension),pphi_i(se_int_n_nodes),sphi_i(se_int_n_nodes))

      do kn_col=1,se_int_n_nodes
        sn_col=element(se_int)%node(kn_col)
        assemble=.false.

        ! =========================================== !
        ! SBIE AND HBIE AT THE SAME COLLOCATION POINT !
        ! =========================================== !

        if (node(sn_col)%dual_is_common) then
          assemble=.true.
          hp=0.
          gp=0.
          mp=0.
          lp=0.
          hm=0.
          gm=0.
          mm=0.
          lm=0.
          x_i=element(se_int)%x_i_hbie(:,kn_col)
          xi_i=element(se_int)%xi_i_hbie(:,kn_col)
          pphi_i=fbem_phi_hybrid(element(se_int)%type_f1,element(se_int)%delta_f,xi_i)
          sphi_i=fbem_phi_hybrid(element(se_int)%type_f2,element(se_int)%delta_f,xi_i)
          do il=1,problem%n
            hp(:,il,il)=hp(:,il,il)+0.5d0*pphi_i
            lp(:,il,il)=lp(:,il,il)-0.5d0*sphi_i
          end do
          if (boundary(sb_int)%class.eq.fbem_boundary_class_cracklike) then
            do il=1,problem%n
              hm(:,il,il)=hm(:,il,il)+0.5d0*pphi_i
              lm(:,il,il)=lm(:,il,il)+0.5d0*sphi_i
            end do
          end if
          if (node(sn_col)%dual.eq.fbem_dual_burton_miller) then
            alpha=node(sn_col)%alpha
            hp=hp+alpha*mp
            gp=gp+alpha*lp
            if (boundary(sb_int)%class.eq.fbem_boundary_class_cracklike) then
              hm=hm+alpha*mm
              gm=gm+alpha*lm
            end if
          end if

        else

          ! ==== !
          ! SBIE !
          ! ==== !

          if ((node(sn_col)%sbie.eq.fbem_sbie).and.(.not.node_freeterm_added(sn_col))) then
            assemble=.true.
            node_freeterm_added(sn_col)=.true.
            hp=0.
            gp=0.
            hm=0.
            gm=0.
            x_i=element(se_int)%x_i_sbie(:,kn_col)
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
                  call fbem_bem_staela2d_sbie_freeterm(n_c_elements,n_set_at_gn         ,t_set_at_gn,geometric_tolerance,nu,c_plus)
                  call fbem_bem_staela2d_sbie_freeterm(n_c_elements,n_set_at_gn_reversed,t_set_at_gn_reversed,geometric_tolerance,nu,c_minus)
                case (3)
                  call fbem_bem_staela3d_sbie_freeterm(n_c_elements,n_set_at_gn         ,t_set_at_gn,geometric_tolerance,nu,c_plus)
                  call fbem_bem_staela3d_sbie_freeterm(n_c_elements,n_set_at_gn_reversed,t_set_at_gn_reversed,geometric_tolerance,nu,c_minus)
              end select
              ! Deallocate temporary variables
              deallocate (n_set_at_gn,t_set_at_gn)
              deallocate (n_set_at_gn_reversed,t_set_at_gn_reversed)
            !
            ! If the node is not in an edge or vertex of the element
            !
            else
              c_plus=0.d0
              do il=1,problem%n
                c_plus(il,il)=0.5d0
              end do
              c_minus=c_plus
            end if
            ! ADD FREE-TERM MATRIX TO INFLUENCE MATRICES
            do il=1,problem%n
              do ik=1,problem%n
                hp(kn_col,il,ik)=hp(kn_col,il,ik)+c_plus(il,ik)
                hm(kn_col,il,ik)=hm(kn_col,il,ik)+c_minus(il,ik)
              end do
            end do

          end if

          ! ======== !
          ! SBIE MCA !
          ! ======== !

          if (node(sn_col)%sbie.eq.fbem_sbie_mca) then
            assemble=.true.
            hp=0.
            gp=0.
            hm=0.
            gm=0.
            x_i=element(se_int)%x_i_sbie_mca(:,kn_col)
            xi_i=element(se_int)%xi_i_sbie_mca(:,kn_col)
            pphi_i=fbem_phi_hybrid(element(se_int)%type_f1,element(se_int)%delta_f,xi_i)
            do il=1,problem%n
              hp(:,il,il)=hp(:,il,il)+0.5d0*pphi_i
              hm(:,il,il)=hm(:,il,il)+0.5d0*pphi_i
            end do
          end if

          ! ====== !
          !  HBIE  !
          ! ====== !

          if (node(sn_col)%hbie.eq.fbem_hbie) then
            assemble=.true.
            mp=0.
            lp=0.
            mm=0.
            lm=0.
            x_i=element(se_int)%x_i_hbie(:,kn_col)
            xi_i=element(se_int)%xi_i_hbie(:,kn_col)
            sphi_i=fbem_phi_hybrid(element(se_int)%type_f2,element(se_int)%delta_f,xi_i)
            do il=1,problem%n
              lp(:,il,il)=lp(:,il,il)-0.5d0*sphi_i
              lm(:,il,il)=lm(:,il,il)+0.5d0*sphi_i
            end do
            ! HBIE
            if (node(sn_col)%dual.eq.0) then
              hp=mp
              gp=lp
              hm=mm
              gm=lm
            end if
            ! BURTON & MILLER FORMULATION
            if (node(sn_col)%dual.eq.fbem_dual_burton_miller) then
              alpha=node(sn_col)%alpha
              hp=hp+alpha*mp
              gp=gp+alpha*lp
              hm=hm+alpha*mm
              gm=gm+alpha*lm
            end if
          end if

        end if

        ! ======== !
        ! ASSEMBLE !
        ! ======== !

        ! The collocation establishes the equations (rows).
        ! The integration establishes the variables (columns).
        if (assemble) then
          if (region(kr)%space.eq.fbem_half_space) then
            if ((x_i(abs(region(kr)%halfspace_n))-region(kr)%halfspace_x).le.(1.d-12/element(se_int)%csize)) then
              select case (region(kr)%halfspace_bc)
                ! u_k=0
                case (0)
                  stop 'u_k=0 half-space not yet'
                ! t_k=0
                case (1)
                  hp=2.*hp
              end select
            end if
          end if
          select case (boundary(sb_int)%coupling)
            case (fbem_boundary_coupling_be)
              select case (boundary(sb_int)%class)
                case (fbem_boundary_class_ordinary)
                  call assemble_bem_staela_equation(kr,sb_int,sb_int_reversion,se_int,se_int_n_nodes,sn_col,1,hp,gp,hm,gm)
                case (fbem_boundary_class_cracklike)
                  call assemble_bem_staela_equation(kr,sb_int,sb_int_reversion,se_int,se_int_n_nodes,sn_col,1,hp,gp,hm,gm)
                  call assemble_bem_staela_equation(kr,sb_int,sb_int_reversion,se_int,se_int_n_nodes,sn_col,2,mp,lp,mm,lm)
              end select
            case (fbem_boundary_coupling_be_fe)
              select case (boundary(sb_int)%class)
                case (fbem_boundary_class_ordinary)
                  call assemble_bem_staela_equation(kr,sb_int,sb_int_reversion,se_int,se_int_n_nodes,sn_col,1,hp,gp,hm,gm)
                case (fbem_boundary_class_cracklike)
                  call assemble_bem_staela_equation(kr,sb_int,sb_int_reversion,se_int,se_int_n_nodes,sn_col,1,hp,gp,hm,gm)
                  call assemble_bem_staela_equation(kr,sb_int,sb_int_reversion,se_int,se_int_n_nodes,sn_col,2,mp,lp,mm,lm)
              end select
            case (fbem_boundary_coupling_be_be,fbem_boundary_coupling_be_fe_be)
              if (sb_int_reversion) then
                call assemble_bem_staela_equation(kr,sb_int,sb_int_reversion,se_int,se_int_n_nodes,sn_col,2,hp,gp,hm,gm)
              else
                call assemble_bem_staela_equation(kr,sb_int,sb_int_reversion,se_int,se_int_n_nodes,sn_col,1,hp,gp,hm,gm)
              end if
          end select
        end if

      end do ! Loop through the NODES of the ELEMENT

      ! Deallocate element-wise data structures
      deallocate (hp,gp,hm,gm)
      deallocate (mp,lp,mm,lm)
      deallocate (xi_i,pphi_i,sphi_i)

    end do ! Loop through the ELEMENTS of the BOUNDARY

  end do ! Loop through the BOUNDARIES of the REGION

  ! ----------------------------------- !
  ! INTEGRATION of REGION BE BODY LOADS !
  ! ----------------------------------- !

  do kb_int=1,region(kr)%n_be_bodyloads
    sb_int=region(kr)%be_bodyload(kb_int)
    sp_int=be_bodyload(sb_int)%part
    if (be_bodyload(sb_int)%coupling.eq.0) cycle
    do ke_int=1,part(sp_int)%n_elements
      se_int=part(sp_int)%element(ke_int)
      se_int_n_nodes=element(se_int)%n_nodes
      se_int_n_snodes=element(se_int)%n_nodes

      allocate (h(se_int_n_nodes,problem%n,problem%n))
      allocate (xi_i(element(se_int)%n_dimension),pphi_i(se_int_n_nodes),sphi_i(se_int_n_nodes))

      do kn_col=1,se_int_n_nodes
        sn_col=element(se_int)%node(kn_col)
        assemble=.false.

        ! ==== !
        ! SBIE !
        ! ==== !

        if ((node(sn_col)%sbie.eq.fbem_sbie).and.(.not.node_freeterm_added(sn_col))) then
          assemble=.true.
          node_freeterm_added(sn_col)=.true.
          h=0.d0
          if (node(sn_col)%sbie_lineload_end_boundary) then
            do il=1,problem%n
              h(kn_col,il,il)=0.5d0
            end do
          else
            do il=1,problem%n
              h(kn_col,il,il)=1.d0
            end do
          end if 
        end if

        ! ======== !
        ! SBIE MCA !
        ! ======== !

        if (node(sn_col)%sbie.eq.fbem_sbie_mca) then
          assemble=.true.
          xi_i=element(se_int)%xi_i_sbie_mca(:,kn_col)
          pphi_i=fbem_phi_hybrid(element(se_int)%type_f1,element(se_int)%delta_f,xi_i)
          h=0.
          do il=1,problem%n
            h(:,il,il)=h(:,il,il)+pphi_i
          end do
        end if

        ! ======== !
        ! ASSEMBLE !
        ! ======== !

        ! The collocation establishes the equations (rows).
        ! The integration establishes the variables (columns).
        if (assemble) then


          !
          ! De momento solo beam line and shell surface coupling
          !
          do il=1,problem%n
            row=node(sn_col)%row(il,1)
            do ik=1,problem%n
              do kn_int=1,se_int_n_nodes
                sn_int=element(se_int)%node(kn_int)
                sn_fe=element(se_int)%element_node(kn_int)
                ! If u_k=u_k_{+}=u_k_{-} is unknown
                if (node(sn_fe)%ctype(ik,1).eq.1) then
                  col=node(sn_fe)%col(ik,1)
                  A_r(row,col)=A_r(row,col)+h(kn_int,il,ik)
                ! If u_k=u_k_{+}=u_k_{-} is known
                else
                  b_r(row,1)=b_r(row,1)-h(kn_int,il,ik)*node(sn_fe)%cvalue_r(ik,1,1)
                end if
              end do
            end do
          end do



        end if

      end do ! Loop through the NODES of the ELEMENT

      ! Deallocate element-wise data structures
      deallocate (h)
      deallocate (xi_i,pphi_i,sphi_i)

    end do ! Loop through the ELEMENTS of the BE BODY LOAD

  end do ! Loop through the BE BODY LOADS of the REGION

  if (verbose_level.ge.2) then
    write(fmtstr,*) '(a24,1x,i',fbem_nchar_int(region(kr)%id),')'
    call fbem_trimall(fmtstr)
    call fbem_timestamp_message(output_unit,2)
    write(output_unit,fmtstr) 'END assembling BE region',region(kr)%id
  end if

end subroutine build_lse_mechanics_bem_staela

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

subroutine build_lse_mechanics_bem_staela_element(kr,sb_int,sb_int_reversion,se_int,se_int_n_nodes,mu,nu)

  ! Fortran 2003 intrinsic module
  use iso_fortran_env

  ! fbem modules
  use fbem_data_structures
  use fbem_string_handling
  use fbem_numerical
  use fbem_symmetry
  use fbem_bem_general
  use fbem_bem_staela2d
  use fbem_bem_staela3d

  ! Module of problem variables
  use problem_variables

  ! No implicit variables
  implicit none

  ! I/O variables
  integer                :: kr
  integer                :: sb_int
  logical                :: sb_int_reversion
  integer                :: se_int
  integer                :: se_int_n_nodes
  real(kind=real64)      :: mu
  real(kind=real64)      :: nu
  ! Local variables
  integer                :: ks
  integer                :: ik
  type(fbem_bem_element) :: se_int_data
  logical                :: se_int_reversion
  integer                :: kb_col, sb_col
  logical                :: sb_col_reversion
  integer                :: sp_col
  integer                :: ke_col, se_col
  integer                :: se_col_n_nodes
  integer                :: kn_col, sn_col
  integer                :: kn, kcp
  real(kind=real64)      :: x_i(problem%n), n_i(problem%n)
  ! Dataset at integration points
  logical                :: node_collocated(n_nodes)
  ! Kernels for SBIE integration
  real(kind=real64)      :: h (se_int_n_nodes,problem%n,problem%n), g (se_int_n_nodes,problem%n,problem%n)
  real(kind=real64)      :: ht(se_int_n_nodes,problem%n,problem%n), gt(se_int_n_nodes,problem%n,problem%n)
  real(kind=real64)      :: hp(se_int_n_nodes,problem%n,problem%n), gp(se_int_n_nodes,problem%n,problem%n)
  real(kind=real64)      :: hm(se_int_n_nodes,problem%n,problem%n), gm(se_int_n_nodes,problem%n,problem%n)
  ! Kernels for HBIE integration
  real(kind=real64)      :: m (se_int_n_nodes,problem%n,problem%n), l (se_int_n_nodes,problem%n,problem%n)
  real(kind=real64)      :: mp(se_int_n_nodes,problem%n,problem%n), lp(se_int_n_nodes,problem%n,problem%n)
  real(kind=real64)      :: mm(se_int_n_nodes,problem%n,problem%n), lm(se_int_n_nodes,problem%n,problem%n)
  ! Multiplier for Dual Burton & Miller formulation
  real(kind=real64)      :: alpha
  ! Associated with symmetry
  real(kind=real64)      :: symconf_m(problem%n), symconf_t(problem%n), symconf_r(problem%n), symconf_s
  logical                :: reversed
  ! Assembling control variable
  logical                :: assemble

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

  !
  ! Loop through SYMMETRICAL ELEMENTS for INTEGRATION
  !
  do ks=1,n_symelements
    ! SYMMETRY SETUP
    call fbem_symmetry_multipliers(ks,problem%n,n_symplanes,symplane_m,symplane_s,symplane_t,symplane_r,&
                                   symconf_m,symconf_s,symconf_t,symconf_r,reversed)
    se_int_reversion=sb_int_reversion.neqv.reversed
    do kn=1,se_int_n_nodes
      se_int_data%x(:,kn)=symconf_m*element(se_int)%x_gn(:,kn)
    end do
    ! INITIALIZE PRECALCULATED DATASETS
    call se_int_data%init_precalculated_datasets(n_precalsets,precalset_gln)

    ! ========================= !
    ! COLLOCATION AT BOUNDARIES !
    ! ========================= !

    !
    ! Loop through the BOUNDARIES of the REGION for COLLOCATION
    !
    node_collocated=.false.
    do kb_col=1,region(kr)%n_boundaries
      sb_col=region(kr)%boundary(kb_col)
      sb_col_reversion=region(kr)%boundary_reversion(kb_col)
      sp_col=boundary(sb_col)%part
      !
      ! Loop through the ELEMENTS of the BOUNDARY for COLLOCATION
      !
      do ke_col=1,part(sp_col)%n_elements
        se_col=part(sp_col)%element(ke_col)
        se_col_n_nodes=element(se_col)%n_nodes
        !
        ! Loop through the NODES of the ELEMENT for COLLOCATION
        !
        do kn_col=1,se_col_n_nodes
          ! COLLOCATION NODE
          sn_col=element(se_col)%node(kn_col)
          assemble=.false.

          ! --------------------------------------------------- !
          ! SBIE & HBIE AT THE SAME NON-NODAL COLLOCATION POINT !
          ! --------------------------------------------------- !

          if (node(sn_col)%dual_is_common) then
            ! INITIALIZE
            assemble=.true.
            x_i=element(se_col)%x_i_hbie(:,kn_col)
            n_i=element(se_col)%n_i_hbie(:,kn_col)
            if (sb_col_reversion) n_i=-n_i
            ! CALCULATE INFLUENCE MATRICES
            select case (problem%n)
              case (2)
                call fbem_bem_staela2d_sbie_auto(se_int_data,se_int_reversion,x_i,mu,nu,qsi_parameters,qsi_ns_max,h,g)
                call fbem_bem_staela2d_hbie_auto(se_int_data,se_int_reversion,x_i,n_i,mu,nu,qsi_parameters,qsi_ns_max,m,l)
              case (3)
                call fbem_bem_staela3d_sbie_auto(se_int_data,se_int_reversion,x_i,mu,nu,qsi_parameters,qsi_ns_max,h,g)
                call fbem_bem_staela3d_hbie_auto(se_int_data,se_int_reversion,x_i,n_i,mu,nu,qsi_parameters,qsi_ns_max,m,l)
            end select
            if (region(kr)%space.eq.fbem_half_space) then
              select case (problem%n)
                case (2)
                  x_i(abs(region(kr)%halfspace_n))=2.d0*region(kr)%halfspace_x-x_i(abs(region(kr)%halfspace_n))
                  n_i(abs(region(kr)%halfspace_n))=-n_i(abs(region(kr)%halfspace_n))
                  call fbem_bem_staela2d_hfc_sbie_auto(se_int_data,se_int_reversion,x_i,mu,nu,region(kr)%halfspace_n,region(kr)%halfspace_x,qsi_parameters,qsi_ns_max,hp,gp)
                  stop 'Error: staela: 2D half-space HBIE not available'
                case (3)
                  call fbem_bem_staela3d_hsc_sbie_auto(se_int_data,se_int_reversion,x_i,mu,nu,region(kr)%halfspace_n,region(kr)%halfspace_x,qsi_parameters,qsi_ns_max,hp,gp)
                  call fbem_bem_staela3d_hsc_hbie_auto(se_int_data,se_int_reversion,x_i,n_i,mu,nu,region(kr)%halfspace_n,region(kr)%halfspace_x,qsi_parameters,qsi_ns_max,mp,lp)
              end select
              select case (region(kr)%halfspace_bc)
                case (0)
                  stop 'Error: staela: region(kr)%halfspace_bc=0 not yet'
                case (1)
                  h=h+hp
                  g=g+gp
                  m=m+mp
                  l=l+lp
              end select
            end if
            ! MODIFY INFLUENCE MATRICES ACCORDING TO SYMMETRY CONFIGURATION
            do ik=1,problem%n
              h(:,:,ik)=symconf_t(ik)*h(:,:,ik)
              g(:,:,ik)=symconf_t(ik)*g(:,:,ik)
              m(:,:,ik)=symconf_t(ik)*m(:,:,ik)
              l(:,:,ik)=symconf_t(ik)*l(:,:,ik)
            end do
            ! BUILD INFLUENCE MATRICES WITH +N and -N
            hp=h
            gp=g
            mp=m
            lp=l
            hm=-h
            gm=g
            mm=-m
            lm=l
            if (node(sn_col)%dual.eq.fbem_dual_burton_miller) then
              alpha=node(sn_col)%alpha
              hp=hp+alpha*mp
              gp=gp+alpha*lp
              hm=hm+alpha*mm
              gm=gm+alpha*lm
            end if

          else

            ! ------------------------------- !
            ! SBIE AT NODAL COLLOCATION POINT !
            ! ------------------------------- !

            if (node(sn_col)%sbie.eq.fbem_sbie.and.(.not.node_collocated(sn_col))) then
              ! INITIALIZE
              assemble=.true.
              node_collocated(sn_col)=.true.
              x_i=element(se_col)%x_i_sbie(:,kn_col)
              ! CALCULATE INFLUENCE MATRICES
              select case (problem%n)
                case (2)
                  call fbem_bem_staela2d_sbie_auto(se_int_data,se_int_reversion,x_i,mu,nu,qsi_parameters,qsi_ns_max,h,g)
                case (3)
                  call fbem_bem_staela3d_sbie_auto(se_int_data,se_int_reversion,x_i,mu,nu,qsi_parameters,qsi_ns_max,h,g)
              end select
              if (region(kr)%space.eq.fbem_half_space) then
                select case (problem%n)
                  case (2)
                    x_i(abs(region(kr)%halfspace_n))=2.d0*region(kr)%halfspace_x-x_i(abs(region(kr)%halfspace_n))
                    n_i(abs(region(kr)%halfspace_n))=-n_i(abs(region(kr)%halfspace_n))
                    call fbem_bem_staela2d_hfc_sbie_auto(se_int_data,se_int_reversion,x_i,mu,nu,region(kr)%halfspace_n,region(kr)%halfspace_x,qsi_parameters,qsi_ns_max,hp,gp)
                  case (3)
                    call fbem_bem_staela3d_hsc_sbie_auto(se_int_data,se_int_reversion,x_i,mu,nu,region(kr)%halfspace_n,region(kr)%halfspace_x,qsi_parameters,qsi_ns_max,hp,gp)
                end select
                select case (region(kr)%halfspace_bc)
                  case (0)
                    stop 'half-space with u_k=0 not available'
                  case (1)
                    h=h+hp
                    g=g+gp
                end select
              end if
              ! MODIFY INFLUENCE MATRICES ACCORDING TO SYMMETRY CONFIGURATION
              do ik=1,problem%n
                h(:,:,ik)=symconf_t(ik)*h(:,:,ik)
                g(:,:,ik)=symconf_t(ik)*g(:,:,ik)
              end do
              ! BUILD INFLUENCE MATRICES WITH +N and -N
              hp=h
              gp=g
              hm=-h
              gm=g
            end if

            ! ----------------------------------- !
            ! SBIE AT NON-NODAL COLLOCATION POINT !
            ! ----------------------------------- !

            if (node(sn_col)%sbie.eq.fbem_sbie_mca) then
              ! INITIALIZE
              assemble=.true.
              x_i=element(se_col)%x_i_sbie_mca(:,kn_col)
              ! CALCULATE INFLUENCE MATRICES
              select case (problem%n)
                case (2)
                  call fbem_bem_staela2d_sbie_auto(se_int_data,se_int_reversion,x_i,mu,nu,qsi_parameters,qsi_ns_max,h,g)
                case (3)
                  call fbem_bem_staela3d_sbie_auto(se_int_data,se_int_reversion,x_i,mu,nu,qsi_parameters,qsi_ns_max,h,g)
              end select
              if (region(kr)%space.eq.fbem_half_space) then
                select case (problem%n)
                  case (2)
                    x_i(abs(region(kr)%halfspace_n))=2.d0*region(kr)%halfspace_x-x_i(abs(region(kr)%halfspace_n))
                    n_i(abs(region(kr)%halfspace_n))=-n_i(abs(region(kr)%halfspace_n))
                    call fbem_bem_staela2d_hfc_sbie_auto(se_int_data,se_int_reversion,x_i,mu,nu,region(kr)%halfspace_n,region(kr)%halfspace_x,qsi_parameters,qsi_ns_max,hp,gp)
                  case (3)
                    call fbem_bem_staela3d_hsc_sbie_auto(se_int_data,se_int_reversion,x_i,mu,nu,region(kr)%halfspace_n,region(kr)%halfspace_x,qsi_parameters,qsi_ns_max,hp,gp)
                end select
                select case (region(kr)%halfspace_bc)
                  case (0)
                    stop 'half-space with u_k=0 not available'
                  case (1)
                    h=h+hp
                    g=g+gp
                end select
              end if
              ! MODIFY INFLUENCE MATRICES ACCORDING TO SYMMETRY CONFIGURATION
              do ik=1,problem%n
                h(:,:,ik)=symconf_t(ik)*h(:,:,ik)
                g(:,:,ik)=symconf_t(ik)*g(:,:,ik)
              end do
              ! BUILD KERNELS WITH N+ AND N-
              hp=h
              gp=g
              hm=-h
              gm= g
            end if

            ! ----------------------------------- !
            ! HBIE AT NON-NODAL COLLOCATION POINT !
            ! ----------------------------------- !

            if (node(sn_col)%hbie.eq.fbem_hbie) then
              ! INITIALIZE
              assemble=.true.
              x_i=element(se_col)%x_i_hbie(:,kn_col)
              n_i=element(se_col)%n_i_hbie(:,kn_col)
              if (sb_col_reversion) n_i=-n_i
              ! CALCULATE INFLUENCE MATRICES
              select case (problem%n)
                case (2)
                  call fbem_bem_staela2d_hbie_auto(se_int_data,se_int_reversion,x_i,n_i,mu,nu,qsi_parameters,qsi_ns_max,m,l)
                case (3)
                  call fbem_bem_staela2d_hbie_auto(se_int_data,se_int_reversion,x_i,n_i,mu,nu,qsi_parameters,qsi_ns_max,m,l)
              end select
              if (region(kr)%space.eq.fbem_half_space) then
                select case (problem%n)
                  case (2)
                    x_i(abs(region(kr)%halfspace_n))=2.d0*region(kr)%halfspace_x-x_i(abs(region(kr)%halfspace_n))
                    n_i(abs(region(kr)%halfspace_n))=-n_i(abs(region(kr)%halfspace_n))
                    stop 'not yet 95'
                  case (3)
                    call fbem_bem_staela3d_hsc_hbie_auto(se_int_data,se_int_reversion,x_i,n_i,mu,nu,region(kr)%halfspace_n,region(kr)%halfspace_x,qsi_parameters,qsi_ns_max,mp,lp)
                end select
                select case (region(kr)%halfspace_bc)
                  case (0)
                    stop 'half-space with u_k=0 not available'
                  case (1)
                    m=m+mp
                    l=l+lp
                end select
              end if
              ! MODIFY INFLUENCE MATRICES ACCORDING TO SYMMETRY CONFIGURATION
              do ik=1,problem%n
                m(:,:,ik)=symconf_t(ik)*m(:,:,ik)
                l(:,:,ik)=symconf_t(ik)*l(:,:,ik)
              end do
              ! BUILD KERNELS WITH N+ AND N-
              mp=m
              lp=l
              mm=-m
              lm=l
              ! HBIE
              if (node(sn_col)%dual.eq.0) then
                hp=mp
                gp=lp
                hm=mm
                gm=lm
              end if
              ! BURTON & MILLER FORMULATION
              if (node(sn_col)%dual.eq.fbem_dual_burton_miller) then
                alpha=node(sn_col)%alpha
                hp=hp+alpha*mp
                gp=gp+alpha*lp
                hm=hm+alpha*mm
                gm=gm+alpha*lm
              end if
            end if

          end if

          ! ======== !
          ! ASSEMBLE !
          ! ======== !

          ! The collocation establishes the equations (rows).
          ! The integration establishes the variables (columns).
          if (assemble) then
            !$omp critical
            select case (boundary(sb_col)%coupling)
              case (fbem_boundary_coupling_be)
                select case (boundary(sb_col)%class)
                  case (fbem_boundary_class_ordinary)
                    call assemble_bem_staela_equation(kr,sb_int,sb_int_reversion,se_int,se_int_n_nodes,sn_col,1,hp,gp,hm,gm)
                  case (fbem_boundary_class_cracklike)
                    call assemble_bem_staela_equation(kr,sb_int,sb_int_reversion,se_int,se_int_n_nodes,sn_col,1,hp,gp,hm,gm)
                    call assemble_bem_staela_equation(kr,sb_int,sb_int_reversion,se_int,se_int_n_nodes,sn_col,2,mp,lp,mm,lm)
                end select
              case (fbem_boundary_coupling_be_fe)
                select case (boundary(sb_col)%class)
                  case (fbem_boundary_class_ordinary)
                    call assemble_bem_staela_equation(kr,sb_int,sb_int_reversion,se_int,se_int_n_nodes,sn_col,1,hp,gp,hm,gm)
                  case (fbem_boundary_class_cracklike)
                    call assemble_bem_staela_equation(kr,sb_int,sb_int_reversion,se_int,se_int_n_nodes,sn_col,1,hp,gp,hm,gm)
                    call assemble_bem_staela_equation(kr,sb_int,sb_int_reversion,se_int,se_int_n_nodes,sn_col,2,mp,lp,mm,lm)
                end select
              case (fbem_boundary_coupling_be_be,fbem_boundary_coupling_be_fe_be)
                if (sb_col_reversion) then
                  call assemble_bem_staela_equation(kr,sb_int,sb_int_reversion,se_int,se_int_n_nodes,sn_col,2,hp,gp,hm,gm)
                else
                  call assemble_bem_staela_equation(kr,sb_int,sb_int_reversion,se_int,se_int_n_nodes,sn_col,1,hp,gp,hm,gm)
                end if
            end select
            !$omp end critical
          end if

        end do ! Loop through the NODES of the ELEMENT for COLLOCATION
      end do ! Loop through the ELEMENTS of the BOUNDARY for COLLOCATION
    end do ! Loop through the BOUNDARIES of the REGION for COLLOCATION

    ! ======================= !
    ! COLLOCATION AT BE LOADS !
    ! ======================= !

    !
    ! Loop through the BE LOADS of the REGION for COLLOCATION
    !
    node_collocated=.false.
    do kb_col=1,region(kr)%n_be_bodyloads
      sb_col=region(kr)%be_bodyload(kb_col)
      sp_col=be_bodyload(sb_col)%part
      select case (be_bodyload(sb_col)%coupling)

        ! ========================
        ! 3D: BEAM - LINE LOAD
        ! ========================

        case (fbem_bl_coupling_beam_line)
          !
          ! Loop through the ELEMENTS of the BE LOAD for COLLOCATION
          !
          do ke_col=1,part(sp_col)%n_elements
            se_col=part(sp_col)%element(ke_col)
            se_col_n_nodes=element(se_col)%n_nodes
            !
            ! Loop through the NODES of the ELEMENT for COLLOCATION
            !
            do kn_col=1,se_col_n_nodes
              ! COLLOCATION NODE
              sn_col=element(se_col)%node(kn_col)
              assemble=.false.

              ! ------------------------------- !
              ! SBIE AT NODAL COLLOCATION POINT !
              ! ------------------------------- !

              if (node(sn_col)%sbie.eq.fbem_sbie.and.(.not.node_collocated(sn_col))) then
                ! INITIALIZE
                assemble=.true.
                node_collocated(sn_col)=.true.
                ! CALCULATE INFLUENCE MATRICES
                h=0
                g=0
                do kcp=1,element(se_col)%cbl_n_cp(kn_col)
                  x_i=element(se_col)%cbl_x_i(:,kcp,kn_col)
                  call fbem_bem_staela3d_sbie_auto(se_int_data,se_int_reversion,x_i,mu,nu,qsi_parameters,qsi_ns_max,ht,gt)
                  if (region(kr)%space.eq.fbem_half_space) then
                    call fbem_bem_staela3d_hsc_sbie_auto(se_int_data,se_int_reversion,x_i,mu,nu,region(kr)%halfspace_n,region(kr)%halfspace_x,qsi_parameters,qsi_ns_max,hp,gp)
                    select case (region(kr)%halfspace_bc)
                      case (0)
                        stop 'half-space with u_k=0 not available'
                      case (1)
                        ht=ht+hp
                        gt=gt+gp
                    end select
                  end if
                  h=h+ht
                  g=g+gt
                end do
                h=h/element(se_col)%cbl_n_cp(kn_col)
                g=g/element(se_col)%cbl_n_cp(kn_col)
                ! MODIFY INFLUENCE MATRICES ACCORDING TO SYMMETRY CONFIGURATION
                do ik=1,problem%n
                  h(:,:,ik)=symconf_t(ik)*h(:,:,ik)
                  g(:,:,ik)=symconf_t(ik)*g(:,:,ik)
                end do
                ! BUILD INFLUENCE MATRICES WITH +N and -N
                hp=h
                gp=g
                hm=-h
                gm=g
              end if

              ! ----------------------------------- !
              ! SBIE AT NON-NODAL COLLOCATION POINT !
              ! ----------------------------------- !

              if (node(sn_col)%sbie.eq.fbem_sbie_mca) then
                ! INITIALIZE
                assemble=.true.
                ! CALCULATE INFLUENCE MATRICES
                h=0
                g=0
                do kcp=1,element(se_col)%cbl_n_cp(kn_col)
                  x_i=element(se_col)%cbl_x_i(:,kcp,kn_col)
                  call fbem_bem_staela3d_sbie_auto(se_int_data,se_int_reversion,x_i,mu,nu,qsi_parameters,qsi_ns_max,ht,gt)
                  if (region(kr)%space.eq.fbem_half_space) then
                    call fbem_bem_staela3d_hsc_sbie_auto(se_int_data,se_int_reversion,x_i,mu,nu,region(kr)%halfspace_n,region(kr)%halfspace_x,qsi_parameters,qsi_ns_max,hp,gp)
                    select case (region(kr)%halfspace_bc)
                      case (0)
                        stop 'half-space with u_k=0 not available'
                      case (1)
                        ht=ht+hp
                        gt=gt+gp
                    end select
                  end if
                  h=h+ht
                  g=g+gt
                end do
                h=h/element(se_col)%cbl_n_cp(kn_col)
                g=g/element(se_col)%cbl_n_cp(kn_col)
                ! MODIFY INFLUENCE MATRICES ACCORDING TO SYMMETRY CONFIGURATION
                do ik=1,problem%n
                  h(:,:,ik)=symconf_t(ik)*h(:,:,ik)
                  g(:,:,ik)=symconf_t(ik)*g(:,:,ik)
                end do
                ! BUILD KERNELS WITH N+ AND N-
                hp=h
                gp=g
                hm=-h
                gm= g
              end if

              ! ======== !
              ! ASSEMBLE !
              ! ======== !

              ! The collocation establishes the equations (rows).
              ! The integration establishes the variables (columns).
              if (assemble) then
                !$omp critical
                select case (be_bodyload(sb_col)%coupling)
                  case (fbem_bl_coupling_beam_line,fbem_bl_coupling_shell_surface)
                    call assemble_bem_staela_equation(kr,sb_int,sb_int_reversion,se_int,se_int_n_nodes,sn_col,1,hp,gp,hm,gm)
                  case default
                    stop 'not yet 97'
                end select
                !$omp end critical
              end if

            end do ! Loop through the NODES of the ELEMENT for COLLOCATION
          end do ! Loop through the ELEMENTS of the BE LOAD for COLLOCATION

        ! ========================
        ! 2D: BEAM - LINE LOAD
        ! 3D: SHELL - SURFACE LOAD
        ! ========================

        case (fbem_bl_coupling_shell_surface)
          !
          ! Loop through the ELEMENTS of the BE LOAD for COLLOCATION
          !
          do ke_col=1,part(sp_col)%n_elements
            se_col=part(sp_col)%element(ke_col)
            se_col_n_nodes=element(se_col)%n_nodes
            !
            ! Loop through the NODES of the ELEMENT for COLLOCATION
            !
            do kn_col=1,se_col_n_nodes
              ! COLLOCATION NODE
              sn_col=element(se_col)%node(kn_col)
              assemble=.false.

              ! ------------------------------- !
              ! SBIE AT NODAL COLLOCATION POINT !
              ! ------------------------------- !

              if (node(sn_col)%sbie.eq.fbem_sbie.and.(.not.node_collocated(sn_col))) then
                ! INITIALIZE
                assemble=.true.
                node_collocated(sn_col)=.true.
                x_i=element(se_col)%x_i_sbie(:,kn_col)
                ! CALCULATE INFLUENCE MATRICES
                select case (problem%n)
                  case (2)
                    call fbem_bem_staela2d_sbie_auto(se_int_data,se_int_reversion,x_i,mu,nu,qsi_parameters,qsi_ns_max,h,g)
                  case (3)
                    call fbem_bem_staela3d_sbie_auto(se_int_data,se_int_reversion,x_i,mu,nu,qsi_parameters,qsi_ns_max,h,g)
                end select
                if (region(kr)%space.eq.fbem_half_space) then
                  select case (problem%n)
                    case (2)
                      x_i(abs(region(kr)%halfspace_n))=2.d0*region(kr)%halfspace_x-x_i(abs(region(kr)%halfspace_n))
                      n_i(abs(region(kr)%halfspace_n))=-n_i(abs(region(kr)%halfspace_n))
                      call fbem_bem_staela2d_hfc_sbie_auto(se_int_data,se_int_reversion,x_i,mu,nu,region(kr)%halfspace_n,region(kr)%halfspace_x,qsi_parameters,qsi_ns_max,hp,gp)
                    case (3)
                      call fbem_bem_staela3d_hsc_sbie_auto(se_int_data,se_int_reversion,x_i,mu,nu,region(kr)%halfspace_n,region(kr)%halfspace_x,qsi_parameters,qsi_ns_max,hp,gp)
                  end select
                  select case (region(kr)%halfspace_bc)
                    case (0)
                      stop 'half-space with u_k=0 not available'
                    case (1)
                      h=h+hp
                      g=g+gp
                  end select
                end if
                ! MODIFY INFLUENCE MATRICES ACCORDING TO SYMMETRY CONFIGURATION
                do ik=1,problem%n
                  h(:,:,ik)=symconf_t(ik)*h(:,:,ik)
                  g(:,:,ik)=symconf_t(ik)*g(:,:,ik)
                end do
                ! BUILD INFLUENCE MATRICES WITH +N and -N
                hp=h
                gp=g
                hm=-h
                gm=g
              end if

              ! ----------------------------------- !
              ! SBIE AT NON-NODAL COLLOCATION POINT !
              ! ----------------------------------- !

              if (node(sn_col)%sbie.eq.fbem_sbie_mca) then
                ! INITIALIZE
                assemble=.true.
                x_i=element(se_col)%x_i_sbie_mca(:,kn_col)
                ! CALCULATE INFLUENCE MATRICES
                select case (problem%n)
                  case (2)
                    call fbem_bem_staela2d_sbie_auto(se_int_data,se_int_reversion,x_i,mu,nu,qsi_parameters,qsi_ns_max,h,g)
                  case (3)
                    call fbem_bem_staela3d_sbie_auto(se_int_data,se_int_reversion,x_i,mu,nu,qsi_parameters,qsi_ns_max,h,g)
                end select
                if (region(kr)%space.eq.fbem_half_space) then
                  select case (problem%n)
                    case (2)
                      x_i(abs(region(kr)%halfspace_n))=2.d0*region(kr)%halfspace_x-x_i(abs(region(kr)%halfspace_n))
                      n_i(abs(region(kr)%halfspace_n))=-n_i(abs(region(kr)%halfspace_n))
                      call fbem_bem_staela2d_hfc_sbie_auto(se_int_data,se_int_reversion,x_i,mu,nu,region(kr)%halfspace_n,region(kr)%halfspace_x,qsi_parameters,qsi_ns_max,hp,gp)
                    case (3)
                      call fbem_bem_staela3d_hsc_sbie_auto(se_int_data,se_int_reversion,x_i,mu,nu,region(kr)%halfspace_n,region(kr)%halfspace_x,qsi_parameters,qsi_ns_max,hp,gp)
                  end select
                  select case (region(kr)%halfspace_bc)
                    case (0)
                      stop 'half-space with u_k=0 not available'
                    case (1)
                      h=h+hp
                      g=g+gp
                  end select
                end if
                ! MODIFY INFLUENCE MATRICES ACCORDING TO SYMMETRY CONFIGURATION
                do ik=1,problem%n
                  h(:,:,ik)=symconf_t(ik)*h(:,:,ik)
                  g(:,:,ik)=symconf_t(ik)*g(:,:,ik)
                end do
                ! BUILD KERNELS WITH N+ AND N-
                hp=h
                gp=g
                hm=-h
                gm= g
              end if

              ! ======== !
              ! ASSEMBLE !
              ! ======== !

              ! The collocation establishes the equations (rows).
              ! The integration establishes the variables (columns).
              if (assemble) then
                !$omp critical
                select case (be_bodyload(sb_col)%coupling)
                  case (fbem_bl_coupling_beam_line,fbem_bl_coupling_shell_surface)
                    call assemble_bem_staela_equation(kr,sb_int,sb_int_reversion,se_int,se_int_n_nodes,sn_col,1,hp,gp,hm,gm)
                  case default
                    stop 'not yet 97'
                end select
                !$omp end critical
              end if

            end do ! Loop through the NODES of the ELEMENT for COLLOCATION
          end do ! Loop through the ELEMENTS of the BE LOAD for COLLOCATION

      end select

    end do ! Loop through the BE LOADS of the REGION for COLLOCATION

  end do ! Loop through SYMMETRICAL ELEMENTS for INTEGRATION

end subroutine build_lse_mechanics_bem_staela_element

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

subroutine build_lse_mechanics_bem_staela_bl(kr,sb_int,se_int,se_int_n_nodes,mu,nu)

  ! Fortran 2003 intrinsic module
  use iso_fortran_env

  ! fbem modules
  use fbem_data_structures
  use fbem_string_handling
  use fbem_numerical
  use fbem_symmetry
  use fbem_bem_general
  use fbem_bem_staela2d
  use fbem_bem_staela3d

  ! Module of problem variables
  use problem_variables

  ! No implicit variables
  implicit none

  ! I/O variables
  integer                :: kr
  integer                :: sb_int
  integer                :: se_int
  integer                :: se_int_n_nodes
  real(kind=real64)      :: mu
  real(kind=real64)      :: nu
  ! Local variables
  integer                :: ks
  integer                :: ik
  type(fbem_bem_element) :: se_int_data
  integer                :: kb_col, sb_col
  logical                :: sb_col_reversion
  integer                :: sp_col
  integer                :: ke_col, se_col
  integer                :: se_col_n_nodes
  integer                :: kn_col, sn_col
  integer                :: kn, kcp
  real(kind=real64), allocatable :: xi_i(:)
  real(kind=real64)      :: x_i(problem%n), n_i(problem%n), ep1(3), ep2(3), ep3(3)
  integer                :: se_fe
  real(kind=real64)      :: A, Iy, Iz, r_integration
  ! Dataset at integration points
  logical                :: node_collocated(n_nodes)
  ! Integrals for SBIE and HBIE integration
  real(kind=real64), allocatable :: g(:,:,:), gt(:,:,:), gp(:,:,:), l(:,:,:)
  ! Multiplier for Dual Burton & Miller formulation
  real(kind=real64)      :: alpha
  ! Associated with symmetry
  real(kind=real64)      :: symconf_m(problem%n), symconf_t(problem%n), symconf_r(problem%n), symconf_s
  logical                :: reversed
  ! Assembling control variable
  logical                :: assemble

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
  allocate (g (se_int_data%n_snodes,problem%n,problem%n))
  allocate (gt(se_int_data%n_snodes,problem%n,problem%n))
  allocate (gp(se_int_data%n_snodes,problem%n,problem%n))
  allocate (l(se_int_data%n_snodes,problem%n,problem%n))

  !
  ! Loop through SYMMETRICAL ELEMENTS for INTEGRATION
  !
  do ks=1,n_symelements
    ! SYMMETRY SETUP
    call fbem_symmetry_multipliers(ks,problem%n,n_symplanes,symplane_m,symplane_s,symplane_t,symplane_r,&
                                   symconf_m,symconf_s,symconf_t,symconf_r,reversed)
    do kn=1,se_int_n_nodes
      se_int_data%x(:,kn)=symconf_m*element(se_int)%x_gn(:,kn)
    end do
    ! INITIALIZE PRECALCULATED DATASETS
    call se_int_data%init_precalculated_datasets(n_precalsets,precalset_gln)

    ! ========================= !
    ! COLLOCATION AT BOUNDARIES !
    ! ========================= !

    !
    ! Loop through the BOUNDARIES of the REGION for COLLOCATION
    !
    node_collocated=.false.
    do kb_col=1,region(kr)%n_boundaries
      sb_col=region(kr)%boundary(kb_col)
      sb_col_reversion=region(kr)%boundary_reversion(kb_col)
      sp_col=boundary(sb_col)%part
      !
      ! Loop through the ELEMENTS of the BOUNDARY for COLLOCATION
      !
      do ke_col=1,part(sp_col)%n_elements
        se_col=part(sp_col)%element(ke_col)
        se_col_n_nodes=element(se_col)%n_nodes
        !
        ! Loop through the NODES of the ELEMENT for COLLOCATION
        !
        do kn_col=1,se_col_n_nodes
          ! COLLOCATION NODE
          sn_col=element(se_col)%node(kn_col)
          assemble=.false.

          ! --------------------------------------------------- !
          ! SBIE & HBIE AT THE SAME NON-NODAL COLLOCATION POINT !
          ! --------------------------------------------------- !

          if (node(sn_col)%dual_is_common) then
            ! INITIALIZE
            assemble=.true.
            x_i=element(se_col)%x_i_hbie(:,kn_col)
            n_i=element(se_col)%n_i_hbie(:,kn_col)
            if (sb_col_reversion) n_i=-n_i
            ! CALCULATE INFLUENCE MATRICES
            select case (problem%n)
              case (2)
                call fbem_bem_staela2d_sbie_bl_auto(se_int_data,x_i,mu,nu,qsi_parameters,qsi_ns_max,g)
                call fbem_bem_staela2d_hbie_bl_auto(se_int_data,x_i,n_i,mu,nu,qsi_parameters,qsi_ns_max,l)
              case (3)
                call fbem_bem_staela3d_sbie_bl_auto(se_int_data,x_i,mu,nu,qsi_parameters,qsi_ns_max,g)
                call fbem_bem_staela3d_hbie_bl_auto(se_int_data,x_i,n_i,mu,nu,qsi_parameters,qsi_ns_max,l)
            end select
            if (region(kr)%space.eq.fbem_half_space) then
              select case (problem%n)
                case (2)
                  x_i(abs(region(kr)%halfspace_n))=2.d0*region(kr)%halfspace_x-x_i(abs(region(kr)%halfspace_n))
                  n_i(abs(region(kr)%halfspace_n))=-n_i(abs(region(kr)%halfspace_n))
                  stop 'Error: staela: 2D half-space HBIE for BE body loads not available'
                case (3)
                  stop 'Error: staela: 3D half-space HBIE for BE body loads not available'
              end select
              select case (region(kr)%halfspace_bc)
                case (0)
                  stop 'Error: staela: region(kr)%halfspace_bc=0 not yet'
                case (1)
              end select
            end if
            ! MODIFY INFLUENCE MATRICES ACCORDING TO SYMMETRY CONFIGURATION
            do ik=1,problem%n
              g(:,:,ik)=symconf_t(ik)*g(:,:,ik)
              l(:,:,ik)=symconf_t(ik)*l(:,:,ik)
            end do
            ! BUILD INFLUENCE MATRICES WITH +N and -N
            if (node(sn_col)%dual.eq.fbem_dual_burton_miller) then
              g=g+node(sn_col)%alpha*l
            end if

          else

            ! ------------------------------- !
            ! SBIE AT NODAL COLLOCATION POINT !
            ! ------------------------------- !

            if (node(sn_col)%sbie.eq.fbem_sbie.and.(.not.node_collocated(sn_col))) then
              ! INITIALIZE
              assemble=.true.
              node_collocated(sn_col)=.true.
              x_i=element(se_col)%x_i_sbie(:,kn_col)
              ! CALCULATE INFLUENCE MATRICES
              select case (problem%n)
                case (2)
                  call fbem_bem_staela2d_sbie_bl_auto(se_int_data,x_i,mu,nu,qsi_parameters,qsi_ns_max,g)
                case (3)
                  call fbem_bem_staela3d_sbie_bl_auto(se_int_data,x_i,mu,nu,qsi_parameters,qsi_ns_max,g)
              end select
              if (region(kr)%space.eq.fbem_half_space) then
                select case (problem%n)
                  case (2)
                    x_i(abs(region(kr)%halfspace_n))=2.d0*region(kr)%halfspace_x-x_i(abs(region(kr)%halfspace_n))
                    n_i(abs(region(kr)%halfspace_n))=-n_i(abs(region(kr)%halfspace_n))
                    stop 'not yet 98'
                  case (3)
                    call fbem_bem_staela3d_hsc_sbie_bl_auto(se_int_data,x_i,mu,nu,region(kr)%halfspace_n,region(kr)%halfspace_x,qsi_parameters,qsi_ns_max,gp)
                end select
                select case (region(kr)%halfspace_bc)
                  case (0)
                    stop 'half-space with u_k=0 not available'
                  case (1)
                    g=g+gp
                end select
              end if
              ! MODIFY INFLUENCE MATRICES ACCORDING TO SYMMETRY CONFIGURATION
              do ik=1,problem%n
                g(:,:,ik)=symconf_t(ik)*g(:,:,ik)
              end do
            end if

            ! ----------------------------------- !
            ! SBIE AT NON-NODAL COLLOCATION POINT !
            ! ----------------------------------- !

            if (node(sn_col)%sbie.eq.fbem_sbie_mca) then
              ! INITIALIZE
              assemble=.true.
              x_i=element(se_col)%x_i_sbie_mca(:,kn_col)
              ! CALCULATE INFLUENCE MATRICES
              select case (problem%n)
                case (2)
                  call fbem_bem_staela2d_sbie_bl_auto(se_int_data,x_i,mu,nu,qsi_parameters,qsi_ns_max,g)
                case (3)
                  call fbem_bem_staela3d_sbie_bl_auto(se_int_data,x_i,mu,nu,qsi_parameters,qsi_ns_max,g)
              end select
              if (region(kr)%space.eq.fbem_half_space) then
                select case (problem%n)
                  case (2)
                    x_i(abs(region(kr)%halfspace_n))=2.d0*region(kr)%halfspace_x-x_i(abs(region(kr)%halfspace_n))
                    n_i(abs(region(kr)%halfspace_n))=-n_i(abs(region(kr)%halfspace_n))
                    stop 'not yet 101'
                  case (3)
                    call fbem_bem_staela3d_hsc_sbie_bl_auto(se_int_data,x_i,mu,nu,region(kr)%halfspace_n,region(kr)%halfspace_x,qsi_parameters,qsi_ns_max,gp)
                end select
                select case (region(kr)%halfspace_bc)
                  case (0)
                    stop 'half-space with u_k=0 not available'
                  case (1)
                    g=g+gp
                end select
              end if
              ! MODIFY INFLUENCE MATRICES ACCORDING TO SYMMETRY CONFIGURATION
              do ik=1,problem%n
                g(:,:,ik)=symconf_t(ik)*g(:,:,ik)
              end do
            end if

            ! ----------------------------------- !
            ! HBIE AT NON-NODAL COLLOCATION POINT !
            ! ----------------------------------- !

            if (node(sn_col)%hbie.eq.fbem_hbie) then
              ! INITIALIZE
              assemble=.true.
              x_i=element(se_col)%x_i_hbie(:,kn_col)
              n_i=element(se_col)%n_i_hbie(:,kn_col)
              if (sb_col_reversion) n_i=-n_i
              ! CALCULATE INFLUENCE MATRICES
              select case (problem%n)
                case (2)
                  call fbem_bem_staela2d_hbie_bl_auto(se_int_data,x_i,n_i,mu,nu,qsi_parameters,qsi_ns_max,l)
                case (3)
                  call fbem_bem_staela3d_hbie_bl_auto(se_int_data,x_i,n_i,mu,nu,qsi_parameters,qsi_ns_max,l)
              end select
              if (region(kr)%space.eq.fbem_half_space) then
                select case (problem%n)
                  case (2)
                    x_i(abs(region(kr)%halfspace_n))=2.d0*region(kr)%halfspace_x-x_i(abs(region(kr)%halfspace_n))
                    n_i(abs(region(kr)%halfspace_n))=-n_i(abs(region(kr)%halfspace_n))
                    stop 'not yet 104'
                  case (3)
                    stop 'not yet 105'
                end select
                select case (region(kr)%halfspace_bc)
                  case (0)
                    stop 'half-space with u_k=0 not available'
                  case (1)
                    stop 'not yet 106'
                end select
              end if
              ! MODIFY INFLUENCE MATRICES ACCORDING TO SYMMETRY CONFIGURATION
              do ik=1,problem%n
                l(:,:,ik)=symconf_t(ik)*l(:,:,ik)
              end do
              ! HBIE
              if (node(sn_col)%dual.eq.0) then
                g=l
              end if
              ! BURTON & MILLER FORMULATION
              if (node(sn_col)%dual.eq.fbem_dual_burton_miller) then
                alpha=node(sn_col)%alpha
                g=g+alpha*l
              end if
            end if

          end if

          ! ======== !
          ! ASSEMBLE !
          ! ======== !

          ! The collocation establishes the equations (rows).
          ! The integration establishes the variables (columns).
          if (assemble) then
            !$omp critical
            select case (boundary(sb_col)%coupling)
              case (fbem_boundary_coupling_be)
                select case (boundary(sb_col)%class)
                  case (fbem_boundary_class_ordinary)
                    call assemble_bem_bl_staela_equation(sb_int,se_int,se_int_data%n_snodes,sn_col,1,g)
                  case (fbem_boundary_class_cracklike)
                    call assemble_bem_bl_staela_equation(sb_int,se_int,se_int_data%n_snodes,sn_col,1,g)
                    call assemble_bem_bl_staela_equation(sb_int,se_int,se_int_data%n_snodes,sn_col,2,l)
                end select
              case (fbem_boundary_coupling_be_fe)
                select case (boundary(sb_col)%class)
                  case (fbem_boundary_class_ordinary)
                    call assemble_bem_bl_staela_equation(sb_int,se_int,se_int_data%n_snodes,sn_col,1,g)
                  case (fbem_boundary_class_cracklike)
                    call assemble_bem_bl_staela_equation(sb_int,se_int,se_int_data%n_snodes,sn_col,1,g)
                    call assemble_bem_bl_staela_equation(sb_int,se_int,se_int_data%n_snodes,sn_col,2,l)
                end select
              case (fbem_boundary_coupling_be_be,fbem_boundary_coupling_be_fe_be)
                if (sb_col_reversion) then
                  call assemble_bem_bl_staela_equation(sb_int,se_int,se_int_data%n_snodes,sn_col,2,g)
                else
                  call assemble_bem_bl_staela_equation(sb_int,se_int,se_int_data%n_snodes,sn_col,1,g)
                end if
            end select
            !$omp end critical
          end if

        end do ! Loop through the NODES of the ELEMENT for COLLOCATION
      end do ! Loop through the ELEMENTS of the BOUNDARY for COLLOCATION
    end do ! Loop through the BOUNDARIES of the REGION for COLLOCATION

    ! ======================= !
    ! COLLOCATION AT BE LOADS !
    ! ======================= !

    !
    ! Loop through the BE LOADS of the REGION for COLLOCATION
    !
    node_collocated=.false.
    do kb_col=1,region(kr)%n_be_bodyloads
      sb_col=region(kr)%be_bodyload(kb_col)
      sp_col=be_bodyload(sb_col)%part
      select case (be_bodyload(sb_col)%coupling)

        ! ========================
        ! 3D: BEAM - LINE LOAD
        ! ========================

        case (fbem_bl_coupling_beam_line)
          !
          ! Loop through the ELEMENTS of the BE LOAD for COLLOCATION
          !
          do ke_col=1,part(sp_col)%n_elements
            se_col=part(sp_col)%element(ke_col)
            se_col_n_nodes=element(se_col)%n_nodes

            !
            ! Loop through the NODES of the ELEMENT for COLLOCATION
            !
            do kn_col=1,se_col_n_nodes
              ! COLLOCATION NODE
              sn_col=element(se_col)%node(kn_col)
              assemble=.false.

              ! ------------------------------- !
              ! SBIE AT NODAL COLLOCATION POINT !
              ! ------------------------------- !

              if (node(sn_col)%sbie.eq.fbem_sbie.and.(.not.node_collocated(sn_col))) then
                ! INITIALIZE
                assemble=.true.
                node_collocated(sn_col)=.true.
                ! CALCULATE INFLUENCE MATRICES
                g=0
                do kcp=1,element(se_col)%cbl_n_cp(kn_col)
                  x_i=element(se_col)%cbl_x_i(:,kcp,kn_col)
                  call fbem_bem_staela3d_sbie_bl_auto(se_int_data,x_i,mu,nu,qsi_parameters,qsi_ns_max,gt)
                  if (region(kr)%space.eq.fbem_half_space) then
                    call fbem_bem_staela3d_hsc_sbie_bl_auto(se_int_data,x_i,mu,nu,region(kr)%halfspace_n,region(kr)%halfspace_x,qsi_parameters,qsi_ns_max,gp)
                    select case (region(kr)%halfspace_bc)
                      case (0)
                        stop 'half-space with u_k=0 not available'
                      case (1)
                        gt=gt+gp
                    end select
                  end if
                  g=g+gt
                end do
                g=g/element(se_col)%cbl_n_cp(kn_col)
                ! MODIFY INFLUENCE MATRICES ACCORDING TO SYMMETRY CONFIGURATION
                do ik=1,problem%n
                  g(:,:,ik)=symconf_t(ik)*g(:,:,ik)
                end do
              end if

              ! ----------------------------------- !
              ! SBIE AT NON-NODAL COLLOCATION POINT !
              ! ----------------------------------- !

              if (node(sn_col)%sbie.eq.fbem_sbie_mca) then
                ! INITIALIZE
                assemble=.true.
                ! CALCULATE INFLUENCE MATRICES
                g=0
                do kcp=1,element(se_col)%cbl_n_cp(kn_col)
                  x_i=element(se_col)%cbl_x_i(:,kcp,kn_col)
                  call fbem_bem_staela3d_sbie_bl_auto(se_int_data,x_i,mu,nu,qsi_parameters,qsi_ns_max,gt)
                  if (region(kr)%space.eq.fbem_half_space) then
                    call fbem_bem_staela3d_hsc_sbie_bl_auto(se_int_data,x_i,mu,nu,region(kr)%halfspace_n,region(kr)%halfspace_x,qsi_parameters,qsi_ns_max,gp)
                    select case (region(kr)%halfspace_bc)
                      case (0)
                        stop 'half-space with u_k=0 not available'
                      case (1)
                        gt=gt+gp
                    end select
                  end if
                  g=g+gt
                end do
                g=g/element(se_col)%cbl_n_cp(kn_col)
                ! MODIFY INFLUENCE MATRICES ACCORDING TO SYMMETRY CONFIGURATION
                do ik=1,problem%n
                  g(:,:,ik)=symconf_t(ik)*g(:,:,ik)
                end do
              end if

              ! ======== !
              ! ASSEMBLE !
              ! ======== !

              ! The collocation establishes the equations (rows).
              ! The integration establishes the variables (columns).
              if (assemble) then
                !$omp critical
                select case (be_bodyload(sb_col)%coupling)
                  case (fbem_bl_coupling_beam_line,fbem_bl_coupling_shell_surface)
                    call assemble_bem_bl_staela_equation(sb_int,se_int,se_int_data%n_snodes,sn_col,1,g)
                  case default
                    stop 'not yet 116'
                end select
                !$omp end critical
              end if

            end do ! Loop through the NODES of the ELEMENT for COLLOCATION
          end do ! Loop through the ELEMENTS of the BE LOAD for COLLOCATION

        ! ========================
        ! 2D: BEAM - LINE LOAD
        ! 3D: SHELL - SURFACE LOAD
        ! ========================

        case (fbem_bl_coupling_shell_surface)
          !
          ! Loop through the ELEMENTS of the BE LOAD for COLLOCATION
          !
          do ke_col=1,part(sp_col)%n_elements
            se_col=part(sp_col)%element(ke_col)
            se_col_n_nodes=element(se_col)%n_nodes
            !
            ! Loop through the NODES of the ELEMENT for COLLOCATION
            !
            do kn_col=1,se_col_n_nodes
              ! COLLOCATION NODE
              sn_col=element(se_col)%node(kn_col)
              assemble=.false.

              ! ------------------------------- !
              ! SBIE AT NODAL COLLOCATION POINT !
              ! ------------------------------- !

              if (node(sn_col)%sbie.eq.fbem_sbie.and.(.not.node_collocated(sn_col))) then
                ! INITIALIZE
                assemble=.true.
                node_collocated(sn_col)=.true.
                x_i=element(se_col)%x_i_sbie(:,kn_col)
                ! CALCULATE INFLUENCE MATRICES
                select case (problem%n)
                  case (2)
                    call fbem_bem_staela2d_sbie_bl_auto(se_int_data,x_i,mu,nu,qsi_parameters,qsi_ns_max,g)
                  case (3)
                    call fbem_bem_staela3d_sbie_bl_auto(se_int_data,x_i,mu,nu,qsi_parameters,qsi_ns_max,g)
                end select
                if (region(kr)%space.eq.fbem_half_space) then
                  select case (problem%n)
                    case (2)
                      x_i(abs(region(kr)%halfspace_n))=2.d0*region(kr)%halfspace_x-x_i(abs(region(kr)%halfspace_n))
                      n_i(abs(region(kr)%halfspace_n))=-n_i(abs(region(kr)%halfspace_n))
                      stop 'not yet 107'
                    case (3)
                      call fbem_bem_staela3d_hsc_sbie_bl_auto(se_int_data,x_i,mu,nu,region(kr)%halfspace_n,region(kr)%halfspace_x,qsi_parameters,qsi_ns_max,gp)
                  end select
                  select case (region(kr)%halfspace_bc)
                    case (0)
                      stop 'half-space with u_k=0 not available'
                    case (1)
                      g=g+gp
                  end select
                end if
                ! MODIFY INFLUENCE MATRICES ACCORDING TO SYMMETRY CONFIGURATION
                do ik=1,problem%n
                  g(:,:,ik)=symconf_t(ik)*g(:,:,ik)
                end do
              end if

              ! ----------------------------------- !
              ! SBIE AT NON-NODAL COLLOCATION POINT !
              ! ----------------------------------- !

              if (node(sn_col)%sbie.eq.fbem_sbie_mca) then
                ! INITIALIZE
                assemble=.true.
                x_i=element(se_col)%x_i_sbie_mca(:,kn_col)
                ! CALCULATE INFLUENCE MATRICES
                select case (problem%n)
                  case (2)
                    call fbem_bem_staela2d_sbie_bl_auto(se_int_data,x_i,mu,nu,qsi_parameters,qsi_ns_max,g)
                  case (3)
                    call fbem_bem_staela3d_sbie_bl_auto(se_int_data,x_i,mu,nu,qsi_parameters,qsi_ns_max,g)
                end select
                if (region(kr)%space.eq.fbem_half_space) then
                  select case (problem%n)
                    case (2)
                      x_i(abs(region(kr)%halfspace_n))=2.d0*region(kr)%halfspace_x-x_i(abs(region(kr)%halfspace_n))
                      n_i(abs(region(kr)%halfspace_n))=-n_i(abs(region(kr)%halfspace_n))
                      stop 'not yet 110'
                    case (3)
                      call fbem_bem_staela3d_hsc_sbie_bl_auto(se_int_data,x_i,mu,nu,region(kr)%halfspace_n,region(kr)%halfspace_x,qsi_parameters,qsi_ns_max,gp)
                  end select
                  select case (region(kr)%halfspace_bc)
                    case (0)
                      stop 'half-space with u_k=0 not available'
                    case (1)
                      g=g+gp
                  end select
                end if
                ! MODIFY INFLUENCE MATRICES ACCORDING TO SYMMETRY CONFIGURATION
                do ik=1,problem%n
                  g(:,:,ik)=symconf_t(ik)*g(:,:,ik)
                end do
              end if

              ! ======== !
              ! ASSEMBLE !
              ! ======== !

              ! The collocation establishes the equations (rows).
              ! The integration establishes the variables (columns).
              if (assemble) then
                !$omp critical
                select case (be_bodyload(sb_col)%coupling)
                  case (fbem_bl_coupling_beam_line,fbem_bl_coupling_shell_surface)
                    call assemble_bem_bl_staela_equation(sb_int,se_int,se_int_data%n_snodes,sn_col,1,g)
                  case default
                    stop 'not yet 116'
                end select
                !$omp end critical
              end if

            end do ! Loop through the NODES of the ELEMENT for COLLOCATION
          end do ! Loop through the ELEMENTS of the BE LOAD for COLLOCATION

      end select ! Switch between BE LOAD COUPLING

    end do ! Loop through the BE LOADS of the REGION for COLLOCATION

  end do ! Loop through SYMMETRICAL ELEMENTS for INTEGRATION

end subroutine build_lse_mechanics_bem_staela_bl
