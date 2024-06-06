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

subroutine build_fem_f_distributed_harmonic(kf,se,ndof_fe,f_e)

  ! Fortran 2003 intrinsic module
  use iso_fortran_env

  ! fbem modules
  use fbem_data_structures
  use fbem_string_handling
  use fbem_numerical
  use fbem_shape_functions
  use fbem_geometry
  use fbem_symmetry
  use fbem_quasisingular_integration
  use fbem_telles_transformation
  use fbem_bem_general
  use fbem_bem_staela2d
  use fbem_bem_staela3d
  use fbem_fem_beams
  use fbem_fem_shells
  use fbem_fem_solids

  ! Module of problem variables
  use problem_variables

  ! No implicit variables
  implicit none

  ! I/O variables
  integer                           :: kf
  integer                           :: se
  integer                           :: ndof_fe
  complex(kind=real64)              :: f_e(ndof_fe)
  ! Local variables
  real(kind=real64)                 :: omega
  integer                           :: ndof_se, ndof_load
  real(kind=real64), allocatable    :: Q(:,:), Ld(:,:), nodal_axes(:,:,:)
  complex(kind=real64), allocatable :: t_e(:)
  integer                           :: sr, kn, kc
  real(kind=real64)                 :: nu
  real(kind=real64)                 :: ep1load(3)

  !
  ! Region of the element (in order to take material properties)
  !
  sr=fe_subregion(part(element(se)%part)%entity)%region

  omega=frequency(kf)

  f_e=0
  select case (element(se)%n_dimension)

    ! ==============================================================================================================================
    ! ONE-DIMENSIONAL ELEMENTS
    ! ==============================================================================================================================

    case (1)

      select case (element(se)%fe_type)

        case (0)

          !-------------------------------------------------------------------------------------------------------------------------
          ! DEGENERATED BEAM FINITE ELEMENT
          !
          if ((element(se)%ctype(1,1).gt.0).or.(element(se)%ctype(2,1).gt.0)) then
            ! Allocate
            ndof_se=3*(problem%n-1)*element(se)%n_nodes
            ndof_load=problem%n*element(se)%n_nodes
            allocate (Q(ndof_se,ndof_load))
            allocate (t_e(ndof_load))
            allocate (Ld(ndof_load,ndof_load))
            ! Calculate
            call fbem_fem_degbeam_Q_midline(problem%n,element(se)%type,element(se)%x_gn,&
                                            element(se)%Q_intmode,element(se)%Q_intngp,Q)
            call fbem_fem_degbeam_L_load(problem%n,element(se)%type,element(se)%v_midnode,Ld)
            ! Build local load if present
            if (element(se)%ctype(1,1).gt.0) then
              t_e=0.d0
              do kn=1,element(se)%n_nodes
                t_e((problem%n*(kn-1)+1):(problem%n*kn))=element(se)%cvalue_c(:,kn,1)
              end do
              f_e(1:ndof_se)=f_e(1:ndof_se)+matmul(Q,matmul(Ld,t_e))
            end if
            ! Build global load if present
            if (element(se)%ctype(2,1).gt.0) then
              t_e=0.d0
              do kn=1,element(se)%n_nodes
                t_e((problem%n*(kn-1)+1):(problem%n*kn))=element(se)%cvalue_c(:,kn,2)
              end do
              f_e(1:ndof_se)=f_e(1:ndof_se)+matmul(Q,t_e)
            end if
            deallocate (Q,Ld,t_e)

!~             allocate (Q(ndof_se,ndof_se))
!~             allocate (t_e(ndof_se))
!~             call fbem_fem_degbeam_Q_body(problem%n,element(se)%type,element(se)%x_gn,element(se)%v_midnode,element(se)%tv_midnode,&
!~                                             element(se)%Q_intmode,element(se)%Q_intngp,Q)
!~             ! Build global load (body) if present
!~             if (element(se)%ctype(2,1).gt.0) then
!~               t_e=0.d0
!~               do kn=1,element(se)%n_nodes
!~                 t_e((3*(problem%n-1)*(kn-1)+1):(3*(problem%n-1)*kn-(2*problem%n-3)))=element(se)%cvalue_c(:,kn,2)
!~               end do
!~               f_e(1:ndof_se)=f_e(1:ndof_se)+matmul(Q,t_e)
!~             end if
!~             deallocate (Q,t_e)

          end if
          !
          !-------------------------------------------------------------------------------------------------------------------------

        case (1,2)

          !-------------------------------------------------------------------------------------------------------------------------
          ! STRAIGHT EULER-BERNOULLI AND TIMOSHENKO BEAM FINITE ELEMENTS
          !
          if ((element(se)%ctype(1,1).gt.0).or.(element(se)%ctype(2,1).gt.0)) then
            !
            ! TO-DO: Los ejes nodales habría que guardarlos para cada nodo y construirlos para cada elemento
            !
            allocate(nodal_axes(problem%n,problem%n,element(se)%n_nodes))
            nodal_axes=0
            do kn=1,element(se)%n_nodes
              do kc=1,problem%n
                nodal_axes(kc,kc,kn)=1
              end do
            end do
            ! Allocate
            ndof_se=3*(problem%n-1)*element(se)%n_nodes
            ndof_load=problem%n*element(se)%n_nodes
            allocate (Q(ndof_se,ndof_load))
            allocate (t_e(ndof_load))
            allocate (Ld(ndof_load,ndof_load))
            !
            ! TODO: Usar el paradigma de region()%material(1)
            ! Use material of the FE REGION
            !
            nu=region(sr)%property_r(3)
            ! Calculate
            call fbem_fem_strbeam_Q_midline(problem%n,element(se)%type,element(se)%fe_options(1),element(se)%fe_type,&
                                                              element(se)%x_gn,element(se)%ep,&
                                                              element(se)%A,element(se)%I(1:3),element(se)%ksh,&
                                                              1.0d0,nu,nodal_axes,Q)
            Ld=fbem_fem_strbeam_L_load(problem%n,element(se)%type,element(se)%ep,nodal_axes)
            ! Build local load if present
            if (element(se)%ctype(1,1).gt.0) then
              t_e=0.d0
              do kn=1,element(se)%n_nodes
                t_e((problem%n*(kn-1)+1):(problem%n*kn))=element(se)%cvalue_c(:,kn,1)
              end do
              f_e(1:ndof_se)=f_e(1:ndof_se)+matmul(Q,matmul(Ld,t_e))
            end if
            ! Build global load if present
            if (element(se)%ctype(2,1).gt.0) then
              t_e=0.d0
              do kn=1,element(se)%n_nodes
                t_e((problem%n*(kn-1)+1):(problem%n*kn))=element(se)%cvalue_c(:,kn,2)
              end do
              f_e(1:ndof_se)=f_e(1:ndof_se)+matmul(Q,t_e)
            end if
            ! Deallocate
            deallocate(Q,Ld,t_e,nodal_axes)
          end if
          !
          !-------------------------------------------------------------------------------------------------------------------------

        case (3)

          !-------------------------------------------------------------------------------------------------------------------------
          ! BAR FINITE ELEMENTS
          !
          if ((element(se)%ctype(1,1).gt.0).or.(element(se)%ctype(2,1).gt.0)) then
            call fbem_error_message(error_unit,0,'element',element(se)%id,'distributed loads not available for bar elements')
          end if
          !
          !-------------------------------------------------------------------------------------------------------------------------

        case (4)

          !-------------------------------------------------------------------------------------------------------------------------
          ! DISCRETE TRANSLATIONAL SPRING FINITE ELEMENTS
          !
          ! DISTRA_EM
          if (element(se)%fe_options(1).eq.1) then
            !
            ! TO-DO: Los ejes nodales habría que guardarlos para cada nodo y construirlos para cada elemento
            !
            allocate(nodal_axes(problem%n,problem%n,2))
            nodal_axes=0
            do kn=1,2
              do kc=1,problem%n
                nodal_axes(kc,kc,kn)=1
              end do
            end do
            !
            allocate (Ld(2*problem%n,2*problem%n),t_e(2*problem%n))
            !
            Ld=fbem_fem_distra_L(problem%n,element(se)%ep,nodal_axes)
            t_e=0
            t_e(1)          =-element(se)%em_Bl*element(se)%em_U/(element(se)%em_R+c_im*omega*element(se)%em_L)
            t_e(problem%n+1)= element(se)%em_Bl*element(se)%em_U/(element(se)%em_R+c_im*omega*element(se)%em_L)
            f_e=matmul(Ld,t_e)
            deallocate(t_e,Ld,nodal_axes)
          end if
          if ((element(se)%ctype(1,1).gt.0).or.(element(se)%ctype(2,1).gt.0)) then
            call fbem_error_message(error_unit,0,'element',element(se)%id,'distributed loads not available for distra elements')
          end if
          !
          !-------------------------------------------------------------------------------------------------------------------------

        case (5)

          !-------------------------------------------------------------------------------------------------------------------------
          ! DISCRETE ROTATIONAL/TRANSLATIONAL SPRING FINITE ELEMENTS
          !
          ! DISROTRA_EM
          if (element(se)%fe_options(1).eq.1) then
            !
            ! TO-DO: Los ejes nodales habría que guardarlos para cada nodo y construirlos para cada elemento
            !
            allocate(nodal_axes(problem%n,problem%n,2))
            nodal_axes=0
            do kn=1,2
              do kc=1,problem%n
                nodal_axes(kc,kc,kn)=1
              end do
            end do
            !
            allocate (Ld(2*3*(problem%n-1),2*3*(problem%n-1)),t_e(2*3*(problem%n-1)))
            !
            Ld=fbem_fem_disrotra_L(problem%n,element(se)%ep,nodal_axes)
            t_e=0
            t_e(1)                =-element(se)%em_Bl*element(se)%em_U/(element(se)%em_R+c_im*omega*element(se)%em_L)
            t_e(3*(problem%n-1)+1)= element(se)%em_Bl*element(se)%em_U/(element(se)%em_R+c_im*omega*element(se)%em_L)
            f_e=matmul(Ld,t_e)
            deallocate(t_e,Ld,nodal_axes)
          end if
          if ((element(se)%ctype(1,1).gt.0).or.(element(se)%ctype(2,1).gt.0)) then
            call fbem_error_message(error_unit,0,'element',element(se)%id,'distributed loads not available for disrotra elements')
          end if
          !
          !-------------------------------------------------------------------------------------------------------------------------

        case (6)

          !-------------------------------------------------------------------------------------------------------------------------
          ! DISCRETE SPRING-DASHPOT
          !
          if ((element(se)%ctype(1,1).gt.0).or.(element(se)%ctype(2,1).gt.0)) then
            call fbem_error_message(error_unit,0,'element',element(se)%id,'distributed loads not available for spring-dashpot elements')
          end if
          !
          !-------------------------------------------------------------------------------------------------------------------------

        case default

          !-------------------------------------------------------------------------------------------------------------------------
          ! OTHER TYPES
          !
          call fbem_error_message(error_unit,0,'element',element(se)%id,'invalid type of 1D element')
          !
          !-------------------------------------------------------------------------------------------------------------------------

      end select

    ! ==============================================================================================================================
    ! TWO-DIMENSIONAL ELEMENTS
    ! ==============================================================================================================================

    case (2)

      select case (problem%n)

        case (2)

          !-------------------------------------------------------------------------------------------------------------------------
          ! SOLID / CONTINUUM ELEMENTS
          !

          !
          ! TODO:Pasar al bucle de ensamblaje de cargas distribuidas
          !
          !
          ! TODO: Falta cargas en este elemento
          !
          !
          !-------------------------------------------------------------------------------------------------------------------------

        case (3)

          !-------------------------------------------------------------------------------------------------------------------------
          ! DEGENERATED SHELL FINITE ELEMENT
          !
          if ((element(se)%ctype(1,1).ne.0).or.(element(se)%ctype(2,1).ne.0)) then
            ! Initialize
            ndof_se=6*element(se)%n_nodes
            ndof_load=3*element(se)%n_nodes
            allocate (Q(ndof_se,ndof_load))
            allocate (Ld(ndof_load,ndof_load))
            allocate (t_e(ndof_load))
            ! Calculate load matrix
            call fbem_fem_degshell_Q_midsurface(element(se)%type,&
                                                element(se)%x_gn,element(se)%v_midnode,element(se)%tv_midnode,&
                                                element(se)%node_n_dof,element(se)%Q_intmode,element(se)%Q_intngp,Q)
            !
            ! TODO: meter t1ref para la carga local
            !
            ep1load=0
            call fbem_fem_degshell_L_load(element(se)%type,element(se)%x_gn,ep1load,geometric_tolerance,Ld)
            ! Load in local coordinates
            if (element(se)%ctype(1,1).eq.1) then
              t_e=0.d0
              do kn=1,element(se)%n_nodes
                t_e((3*(kn-1)+1):(3*kn))=element(se)%cvalue_c(:,kn,1)
              end do
              f_e(1:ndof_se)=f_e(1:ndof_se)+matmul(Q,matmul(Ld,t_e))
            end if
            ! Load in global coordinates
            if (element(se)%ctype(2,1).eq.1) then
              t_e=0.d0
              do kn=1,element(se)%n_nodes
                t_e((3*(kn-1)+1):(3*kn))=element(se)%cvalue_c(:,kn,2)
              end do
              f_e(1:ndof_se)=f_e(1:ndof_se)+matmul(Q,t_e)
            end if
            ! Deallocate
            deallocate (t_e,Ld,Q)
          end if
          !
          !-------------------------------------------------------------------------------------------------------------------------

      end select

    ! ==============================================================================================================================
    ! THREE-DIMENSIONAL ELEMENTS
    ! ==============================================================================================================================

    case (3)

      !-----------------------------------------------------------------------------------------------------------------------------
      ! SOLID / CONTINUUM ELEMENTS
      !
      call fbem_error_message(error_unit,0,'element',element(se)%id,'3D elements not available yet')
      !
      !-----------------------------------------------------------------------------------------------------------------------------

  end select

end subroutine build_fem_f_distributed_harmonic
