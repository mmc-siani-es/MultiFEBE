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

subroutine build_fem_K_static(se,ndof_Kinout,Kinout)

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
  use fbem_fem_beams
  use fbem_fem_shells
  use fbem_fem_solids

  ! Module of problem variables
  use problem_variables

  ! No implicit variables
  implicit none

  ! I/O variables
  integer                        :: se
  integer                        :: ndof_Kinout
  real(kind=real64)              :: Kinout(ndof_Kinout,ndof_Kinout)
  ! Local variables
  integer                        :: sr, kn, kc
  real(kind=real64)              :: E, nu

  real(kind=real64)              :: E11, E22, E33, G12, G13, G23       ! Properties of a orthotropic material
  real(kind=real64)              :: nu12, nu21, nu13, nu31, nu23, nu32 ! Properties of a orthotropic material
  real(kind=real64), allocatable :: nodal_axes(:,:,:)
  integer                        :: ndof_K
  real(kind=real64), allocatable :: K(:,:)

  !
  ! Region of the element (in order to take material properties)
  !
  sr=fe_subregion(part(element(se)%part)%entity)%region
  !
  ! TODO: Usar el paradigma de region()%material(1)
  ! Use material of the FE REGION
  !
  if (region(sr)%subtype.eq.fbem_elastic_orthotropic) then
    ! ORTHOTROPIC
    stop 'orthotropic material not yet available'
    E11 =region(sr)%property_r( 1)
    E22 =region(sr)%property_r( 2)
    E33 =region(sr)%property_r( 3)
    G12 =region(sr)%property_r( 4)
    G13 =region(sr)%property_r( 5)
    G23 =region(sr)%property_r( 6)
    nu12=region(sr)%property_r( 7)
    nu21=region(sr)%property_r( 8)
    nu13=region(sr)%property_r( 9)
    nu31=region(sr)%property_r(10)
    nu23=region(sr)%property_r(11)
    nu32=region(sr)%property_r(12)
  else
    ! ISOTROPIC
    E =region(sr)%property_r(5)
    nu=region(sr)%property_r(3)
  end if

  Kinout=0
  select case (element(se)%n_dimension)

    ! ====================================================================================================================
    ! ONE-DIMENSIONAL ELEMENTS
    ! ====================================================================================================================

    case (1)

      select case (element(se)%fe_type)

        case (0)

          !---------------------------------------------------------------------------------------------------------------
          ! DEGENERATED BEAM FINITE ELEMENT
          !
          ndof_K=3*(problem%n-1)*element(se)%n_nodes
          allocate (K(ndof_K,ndof_K))
          K=0
          select case (problem%n)
            case (2)
              select case (problem%subtype)
                case (fbem_mechanics_plane_strain)
                  call fbem_fem_degbeam_K_static(2,1,&
                                                 element(se)%type,&
                                                 element(se)%x_gn,&
                                                 element(se)%v_midnode,&
                                                 element(se)%tv_midnode,&
                                                 E,nu,element(se)%ksh,element(se)%K_intmode,element(se)%K_intngp,K)
                case (fbem_mechanics_plane_stress)
                  call fbem_fem_degbeam_K_static(2,2,&
                                                 element(se)%type,&
                                                 element(se)%x_gn,&
                                                 element(se)%v_midnode,&
                                                 element(se)%tv_midnode,&
                                                 E,nu,element(se)%ksh,element(se)%K_intmode,element(se)%K_intngp,K)
                case default
                  call fbem_error_message(error_unit,0,'element',element(se)%id,'invalid type of 2D analysis')
              end select
            case (3)
              call fbem_fem_degbeam_K_static(3,0,&
                                             element(se)%type,&
                                             element(se)%x_gn,&
                                             element(se)%v_midnode,&
                                             element(se)%tv_midnode,&
                                             E,nu,element(se)%ksh,element(se)%K_intmode,element(se)%K_intngp,K)
            case default
              call fbem_error_message(error_unit,0,'element',element(se)%id,'invalid problem%n')
          end select
          !
          !---------------------------------------------------------------------------------------------------------------

        case (1,2)

          !---------------------------------------------------------------------------------------------------------------
          ! STRAIGHT EULER-BERNOULLI AND TIMOSHENKO BEAM FINITE ELEMENTS
          !
          ndof_K=3*(problem%n-1)*element(se)%n_nodes
          allocate (K(ndof_K,ndof_K))
          K=0
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
          !
          ! NOTA: viga plane strain (A debe introducirse ya con espesor unitario, y para plane stress con el espesor que toque)
          call fbem_fem_strbeam_K_static(problem%n,&
                                        element(se)%type,element(se)%fe_options(1),element(se)%fe_type,&
                                        element(se)%x_gn,&
                                        element(se)%ep,element(se)%A,element(se)%I(1:3),element(se)%ksh,&
                                        E,nu,nodal_axes,K)
          !
          !---------------------------------------------------------------------------------------------------------------

        case (3)

          !---------------------------------------------------------------------------------------------------------------
          ! BAR FINITE ELEMENTS
          !
          ndof_K=problem%n*element(se)%n_nodes
          allocate (K(ndof_K,ndof_K))
          K=0
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
          !
          call fbem_fem_bar_K_static(problem%n,element(se)%x_gn,element(se)%A,nodal_axes,E,K)
          !
          !---------------------------------------------------------------------------------------------------------------

        case (4)

          !---------------------------------------------------------------------------------------------------------------
          ! DISCRETE TRANSLATIONAL SPRING FINITE ELEMENTS
          !
          ndof_K=2*problem%n
          allocate (K(ndof_K,ndof_K))
          K=0
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
          call fbem_fem_distra_K_static(problem%n,element(se)%ep,nodal_axes,element(se)%k_r,K)
          !
          !---------------------------------------------------------------------------------------------------------------

        case (5)

          !---------------------------------------------------------------------------------------------------------------
          ! DISCRETE ROTATIONAL/TRANSLATIONAL SPRING FINITE ELEMENTS
          !
          ndof_K=3*(problem%n-1)*2
          allocate (K(ndof_K,ndof_K))
          K=0
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
          call fbem_fem_disrotra_K_static(problem%n,element(se)%ep,nodal_axes,element(se)%k_r,K)
          !
          !---------------------------------------------------------------------------------------------------------------

        case default

          !---------------------------------------------------------------------------------------------------------------
          ! OTHER TYPES
          !
          call fbem_error_message(error_unit,0,'element',element(se)%id,'invalid type of 1D element')
          !
          !---------------------------------------------------------------------------------------------------------------

      end select

    ! ====================================================================================================================
    ! TWO-DIMENSIONAL ELEMENTS
    ! ====================================================================================================================

    case (2)

      select case (problem%n)

        case (2)

          !---------------------------------------------------------------------------------------------------------------
          ! SOLID / CONTINUUM ELEMENTS
          !
          ndof_K=2*element(se)%n_nodes
          allocate (K(ndof_K,ndof_K))
          K=0
          select case (problem%subtype)
            case (fbem_mechanics_plane_strain)
              call fbem_fem_solid2d_K_static(element(se)%type,1,&
                                             element(se)%x_gn,element(se)%tn_midnode(3,:),&
                                             E,nu,element(se)%K_intmode,element(se)%K_intngp,K)
            case (fbem_mechanics_plane_stress)
              call fbem_fem_solid2d_K_static(element(se)%type,2,&
                                             element(se)%x_gn,element(se)%tn_midnode(3,:),&
                                             E,nu,element(se)%K_intmode,element(se)%K_intngp,K)
            case default
              call fbem_error_message(error_unit,0,'element',element(se)%id,'invalid type of 2D analysis')
          end select
          !
          !---------------------------------------------------------------------------------------------------------------

        case (3)

          !---------------------------------------------------------------------------------------------------------------
          ! DEGENERATED SHELL FINITE ELEMENT
          !
          ndof_K=6*element(se)%n_nodes
          allocate (K(ndof_K,ndof_K))
          K=0
          call fbem_fem_degshell_K_static(element(se)%type,element(se)%mitc,&
                                          element(se)%x_gn,element(se)%v_midnode,element(se)%tv_midnode,&
                                          element(se)%node_n_dof,E,nu,element(se)%ksh,&
                                          element(se)%K_intmode,element(se)%K_intngp,K)
          !
          !---------------------------------------------------------------------------------------------------------------

      end select

    ! ====================================================================================================================
    ! THREE-DIMENSIONAL ELEMENTS
    ! ====================================================================================================================

    case (3)

      !-------------------------------------------------------------------------------------------------------------------
      ! SOLID / CONTINUUM ELEMENTS
      !
      call fbem_error_message(error_unit,0,'element',element(se)%id,'3D elements not available yet')
      !ndof_K=3*element(se)%n_nodes
      !allocate (K(ndof_K,ndof_K))
      !K=0
      !call fbem_fem_solid3d_K_static(element(se)%type,element(se)%x_gn,E,nu,element(se)%K_intmode,element(se)%K_intngp,K)
      !
      !-------------------------------------------------------------------------------------------------------------------

  end select

  Kinout(1:ndof_K,1:ndof_K)=K

  deallocate(K)

end subroutine build_fem_K_static
