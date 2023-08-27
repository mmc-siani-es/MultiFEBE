! ---------------------------------------------------------------------
! Copyright (C) 2014-2023 Universidad de Las Palmas de Gran Canaria:
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

subroutine build_fem_N(se,se_dim,xi,ndof_u,ndof_Ninout,Ninout)

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
  integer                        :: se_dim
  real(kind=real64)              :: xi(se_dim)
  integer                        :: ndof_u
  integer                        :: ndof_Ninout
  real(kind=real64)              :: Ninout(ndof_u,ndof_Ninout)
  ! Local variables
  integer                        :: sr, kn, kc
  real(kind=real64)              :: xi3d(3), xi1d
  real(kind=real64)              :: E, nu
  real(kind=real64), allocatable :: nodal_axes(:,:,:)
  integer                        :: ndof_N
  real(kind=real64), allocatable :: Lc(:,:), N(:,:), N5(:,:)


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
  else
    ! ISOTROPIC
    E =region(sr)%property_r(5)
    nu=region(sr)%property_r(3)
  end if

  Ninout=0
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
          stop 'not yet 1 build_fem_N.f90'
          !
          !---------------------------------------------------------------------------------------------------------------

        case (1,2)

          !---------------------------------------------------------------------------------------------------------------
          ! STRAIGHT EULER-BERNOULLI AND TIMOSHENKO BEAM FINITE ELEMENTS
          !
          ndof_N=3*(problem%n-1)*element(se)%n_nodes
          allocate (N(ndof_N,ndof_N))
          xi1d=xi(1)
          N=0
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
          N=fbem_fem_strbeam_N(problem%n,element(se)%type,element(se)%fe_options(1),element(se)%fe_type,&
                                  element(se)%x_gn,element(se)%ep,element(se)%A,element(se)%I(1:3),element(se)%ksh,&
                                  E,nu,nodal_axes,xi1d)
          !
          !---------------------------------------------------------------------------------------------------------------

        case (3)

          !---------------------------------------------------------------------------------------------------------------
          ! BAR FINITE ELEMENTS
          !
          stop 'not yet 2 build_fem_N.f90'
          !
          !---------------------------------------------------------------------------------------------------------------

        case (4)

          !---------------------------------------------------------------------------------------------------------------
          ! DISCRETE TRANSLATIONAL SPRING FINITE ELEMENTS
          !
          stop 'not yet 3 build_fem_N.f90'
          !
          !---------------------------------------------------------------------------------------------------------------

        case (5)

          !---------------------------------------------------------------------------------------------------------------
          ! DISCRETE ROTATIONAL/TRANSLATIONAL SPRING FINITE ELEMENTS
          !
          stop 'not yet 4 build_fem_N.f90'
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
          stop 'not yet 5 build_fem_N.f90'
          !
          !---------------------------------------------------------------------------------------------------------------

        case (3)

          !---------------------------------------------------------------------------------------------------------------
          ! DEGENERATED SHELL FINITE ELEMENT
          !
          ndof_N=5*element(se)%n_nodes
          allocate (N5(3,ndof_N))
          N5=0
          xi3d(1)=xi(1)
          xi3d(2)=xi(2)
          xi3d(3)=0
          call fbem_fem_degshell_N(element(se)%type,element(se)%v_midnode,element(se)%tv_midnode,xi3d,N5)
          ! Calculate ENLARGED COORDINATE TRANSFORMATION MATRIX (Lc)
          ndof_N=sum(element(se)%node_n_dof)
          allocate (Lc(ndof_N,ndof_N),N(3,ndof_N))
          call fbem_fem_degshell_Lc_enlarged(element(se)%type,ndof_N,element(se)%node_n_dof,element(se)%tv_midnode,Lc)
          N=matmul(N5,transpose(Lc))
          deallocate(Lc,N5)
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
      stop 'not yet 6 build_fem_N.f90'!
      !-------------------------------------------------------------------------------------------------------------------

  end select

  Ninout(1:ndof_u,1:ndof_N)=N

  deallocate(N)

end subroutine build_fem_N
