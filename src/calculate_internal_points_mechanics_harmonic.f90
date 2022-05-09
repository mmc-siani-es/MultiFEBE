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


subroutine calculate_internal_points_mechanics_harmonic(kf)

  ! Fortran 2003 intrinsic module
  use iso_fortran_env

  ! fbem modules
  use fbem_data_structures
  use fbem_string_handling
  use fbem_shape_functions

  ! Module of problem variables
  use problem_variables

  ! No implicit variables
  implicit none

  ! I/O variables
  integer                           :: kf
  ! Local variables
  integer                           :: kr, kip, kp, ke, se, k_start, k_end, etype, i, j
  complex(kind=real64), allocatable :: value_c(:,:,:)
  real(kind=real64)                 :: aux(10), xi1d, delta_f
  real(kind=real64), allocatable    :: xi(:), phi(:)

  ! Message
  if (verbose_level.ge.1)  write(output_unit,'(a)') 'Calculating internal points solutions at BE regions ...'

  ! Initialize internal points
  do kip=1,n_internalpoints
    internalpoint(kip)%value_c=(0.d0,0.d0)
  end do

  ! Calculate internal points
  do kr=1,n_regions
    if ((region(kr)%class.eq.fbem_be).and.(region(kr)%n_internalpoints.gt.0)) then
      select case (region(kr)%type)
        case (fbem_potential)
          call calculate_internal_points_mechanics_bem_harpot(kf,kr)
        case (fbem_viscoelastic)
          call calculate_internal_points_mechanics_bem_harela(kf,kr)
        case (fbem_poroelastic)
          call calculate_internal_points_mechanics_bem_harpor(kf,kr)
      end select
    end if
  end do

  ! Ending message
  if (verbose_level.ge.1) write(output_unit,'(a)') 'done.'






  if (internalelements) then

    ! Message
    if (verbose_level.ge.1)  write(output_unit,'(a)') 'Map the internal points corresponding to internal elements ...'


    !
    ! Falta el campo incidente...
    !


    do kp=1,internalelements_mesh%n_parts
      kr=internalelements_mesh%part(kp)%entity
      if (kr.eq.0) cycle
      if (region(kr)%class.eq.fbem_be) then
        select case (region(kr)%type)
          case (fbem_potential)
              k_start=1
              k_end  =1
          case (fbem_viscoelastic)
              k_start=1
              k_end  =problem%n
          case (fbem_poroelastic)
              k_start=0
              k_end  =problem%n
        end select
        do ke=1,internalelements_mesh%part(kp)%n_elements
          se=internalelements_mesh%part(kp)%element(ke)
          etype=internalelements_mesh%element(se)%type
          delta_f=internalelements_mesh%element(se)%delta_f
          internalelements_mesh%element(se)%value_c=0
          allocate (xi(internalelements_mesh%element(se)%n_dimension))
          allocate (phi(internalelements_mesh%element(se)%n_nodes))
          allocate (value_c(k_start:k_end,internalelements_mesh%element(se)%n_nodes,0:problem%n))
          ! Copy the corresponding internal point solutions to value_c
          value_c=0
          do j=1,internalelements_mesh%element(se)%n_nodes
            kip=internalelements_mesh%element(se)%internalpoint(j)
            value_c(:,j,:)=internalpoint(kip)%value_c
          end do
          ! Calculate values at element nodes
          select case (internalelements_mesh%element(se)%n_dimension)
            case (1)
              do j=1,internalelements_mesh%element(se)%n_nodes
#               define node j
#               define delta 0.d0
#               define xi xi1d
#               include <xi_1d_at_node.rc>
#               undef delta
#               undef node
#               define delta delta_f
#               include <phi_1d.rc>
#               undef xi
#               undef delta
                do i=1,internalelements_mesh%element(se)%n_nodes
                  internalelements_mesh%element(se)%value_c(:,j,:)=internalelements_mesh%element(se)%value_c(:,j,:)+phi(i)*value_c(:,j,:)
                end do
              end do
            case (2)
              do j=1,internalelements_mesh%element(se)%n_nodes
#               define node j
#               define delta 0.d0
#               include <xi_2d_at_node.rc>
#               undef delta
#               undef node
#               define delta delta_f
#               include <phi_2d.rc>
#               undef delta
                do i=1,internalelements_mesh%element(se)%n_nodes
                  internalelements_mesh%element(se)%value_c(:,j,:)=internalelements_mesh%element(se)%value_c(:,j,:)+phi(i)*value_c(:,j,:)
                end do
              end do
            case (3)
              do j=1,internalelements_mesh%element(se)%n_nodes
#               define node j
#               define delta 0.d0
#               include <xi_3d_at_node.rc>
#               undef delta
#               undef node
#               define delta delta_f
#               include <phi_3d.rc>
#               undef delta
                do i=1,internalelements_mesh%element(se)%n_nodes
                  internalelements_mesh%element(se)%value_c(:,j,:)=internalelements_mesh%element(se)%value_c(:,j,:)+phi(i)*value_c(:,j,:)
                end do
              end do
          end select
          deallocate(xi,phi,value_c)
        end do
      end if
    end do






    ! Ending message
    if (verbose_level.ge.1) write(output_unit,'(a)') 'done.'

  end if







end subroutine calculate_internal_points_mechanics_harmonic
