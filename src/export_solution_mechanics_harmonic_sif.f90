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

subroutine export_solution_mechanics_harmonic_sif(kf,output_fileunit)

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

  ! Frequency
  omega=frequency(kf)
  if (frequency_units.eq.'f') omega=omega*c_1_2pi


  ! ------
  ! Export
  ! ------

  ! Output format
  write(fmtstr,*) '(',fmt_integer,',',2*problem%n+1,fmt_real,')'
  call fbem_trim2b(fmtstr)

  ! Loop through nodes where the SIF must be calculated
  do kn=1,n_nodes
    if (node(kn)%export_sif) then

      ! Only in 2D
      if (problem%n.eq.3) stop 'not implemented yet'

      ! Element of the node
      ke=node(kn)%element(1)
      select case (element(ke)%type)
        !
        ! LINE3
        !
        case (fbem_line3)
          ! Local index of the node in its element
          knl=node(kn)%element_node_iid(1)
          ! Only for quadratic elements
          midnode=element(ke)%node(3)
          select case (knl)
            case (1)
              endnode=element(ke)%node(2)
            case (2)
              endnode=element(ke)%node(1)
          end select
          ! Crack tip local system of coordinates
          do k=1,problem%n
            n_tip(k)=element(ke)%n_gn(k,knl)
            t1_tip(k)=element(ke)%t1_gn(k,knl)
            if (problem%n.eq.3) t2_tip(k)=element(ke)%t2_gn(k,knl)
          end do
          ! Calculate lengths between the tip and the midnode, and the end node.
          select case (knl)
            case (1)
              r_midnode=fbem_length2d_subdivision(element(ke)%type,element(ke)%x_gn,-1.0d0,0.0d0,0,0.0d0,0)
            case (2)
              r_midnode=fbem_length2d_subdivision(element(ke)%type,element(ke)%x_gn, 0.0d0,1.0d0,0,0.0d0,0)
            case default
              stop 'crack-tip element not valid'
          end select
          r_endnode=fbem_length2d_subdivision(element(ke)%type,element(ke)%x_gn,-1.0d0,1.0d0,0,0.0d0,0)
          ratio=r_midnode/r_endnode
          ! Boundary of the node/element
          kb=part(node(kn)%part(1))%entity
          ! Region of the boundary
          kr=boundary(kb)%region(1)
          !
          ! Switch depending on the type of region of the node
          !
          select case (region(kr)%type)

            ! -------------- !
            ! INVISCID FLUID !
            ! -------------- !

            case (fbem_potential)
              stop 'not implemented yet'

            ! ------------------ !
            ! VISCOELASTIC SOLID !
            ! ------------------ !

            case (fbem_viscoelastic)
              select case (region(kr)%subtype)
                !
                ! Viscoelastic or elastic
                ! M.H. Aliabadi, The Boundary Element Method: Volume 2 Applications in solids and structures, 2002. Pages 346-?
                ! J. Dominguez, Boundary Elements in Dynamics, 1993. Pages 509-?
                !
                case (0,fbem_elastic)
                  ! Switch depending on the boundary class
                  select case (boundary(kb)%class)
                    !
                    ! Crack-like boundary
                    !
                    case (fbem_boundary_class_cracklike)
                      mu=region(kr)%property_c(2)
                      kappa=3.d0-4.d0*region(kr)%property_r(3)
                      !kappa=(3.d0-region(kr)%property_r(3))/(1.d0+region(kr)%property_r(3))
                      ! Calculate delta u
                      delta_u_I_midnode=0.d0
                      delta_u_II_midnode=0.d0
                      delta_u_I_endnode=0.d0
                      delta_u_II_endnode=0.d0
                      do k=1,problem%n
                        delta_u_I_midnode =delta_u_I_midnode +(node(midnode)%value_c(k,2)-node(midnode)%value_c(k,1))*n_tip(k)
                        delta_u_II_midnode=delta_u_II_midnode+(node(midnode)%value_c(k,2)-node(midnode)%value_c(k,1))*t1_tip(k)
                        delta_u_I_endnode =delta_u_I_endnode +(node(endnode)%value_c(k,2)-node(endnode)%value_c(k,1))*n_tip(k)
                        delta_u_II_endnode=delta_u_II_endnode+(node(endnode)%value_c(k,2)-node(endnode)%value_c(k,1))*t1_tip(k)
                      end do
                      ! Calculate K_I and K_II using extrapolation of K_I and K_II at mid node and end node.
                      !K_I =mu/(kappa+1.0d0)*dsqrt(c_2pi)/(ratio-1.d0)*(ratio*(delta_u_I_endnode )/dsqrt(r_endnode)-(delta_u_I_midnode )/dsqrt(r_midnode))
                      !K_II=mu/(kappa+1.0d0)*dsqrt(c_2pi)/(ratio-1.d0)*(ratio*(delta_u_II_endnode)/dsqrt(r_endnode)-(delta_u_II_midnode)/dsqrt(r_midnode))
                      ! Calculate K_I and K_II using the end node solution
                      !K_I =mu/(kappa+1.0d0)*dsqrt(c_2pi/r_endnode)*delta_u_I_endnode
                      !K_II=mu/(kappa+1.0d0)*dsqrt(c_2pi/r_endnode)*delta_u_II_endnode
                      ! Calculate K_I and K_II using the mid node solution
                      K_I =mu/(kappa+1.0d0)*dsqrt(c_2pi)*delta_u_I_midnode/dsqrt(r_midnode)
                      K_II=mu/(kappa+1.0d0)*dsqrt(c_2pi)*delta_u_II_midnode/dsqrt(r_midnode)
                      ! Write results
                      write(output_fileunit,fmtstr) node(kn)%id, omega, zabs(K_I), zabs(K_II), fbem_zarg(K_I), fbem_zarg(K_II)
                    !
                    ! Ordinary boundary
                    !
                    case (fbem_boundary_class_ordinary)
                      ! Calculate t tangential and normal
                      t_I=0.d0
                      t_II=0.d0
                      do k=1,problem%n
                        t_I =t_I +node(kn)%value_c(problem%n+k,1)*n_tip(k)
                        t_II=t_II+node(kn)%value_c(problem%n+k,1)*t1_tip(k)
                      end do
                      ! Calculate K_I and K_II using the mid node solution
                      K_I =sqrt(c_2pi*r_endnode)*t_I
                      K_II=sqrt(c_2pi*r_endnode)*t_II
                      ! Write results
                      write(output_fileunit,fmtstr) node(kn)%id, omega, zabs(K_I), zabs(K_II), fbem_zarg(K_I), fbem_zarg(K_II)
                  end select

                !
                ! Other cases (Bardet, etc..)
                !
                case default
                 stop 'not implemented yet'
              end select

            ! ----------------- !
            ! POROELASTIC MEDIUM !
            ! ----------------- !

            case (fbem_poroelastic)
                  ! Switch depending on the boundary class
                  select case (boundary(kb)%class)

                    !
                    ! Crack-like boundary
                    !
                    case (fbem_boundary_class_cracklike)
                      mu=region(kr)%property_c(4)
                      kappa=3.0d0-4.0d0*region(kr)%property_r(6)
                      !kappa=(3.d0-region(kr)%property_r(6))/(1.d0+region(kr)%property_r(6))
                      ! Calculate delta u
                      delta_u_I_midnode=0.d0
                      delta_u_II_midnode=0.d0
                      delta_u_I_endnode=0.d0
                      delta_u_II_endnode=0.d0
                      do k=1,problem%n
                        delta_u_I_midnode =delta_u_I_midnode +(node(midnode)%value_c(k,2)-node(midnode)%value_c(k,1))*n_tip(k)
                        delta_u_II_midnode=delta_u_II_midnode+(node(midnode)%value_c(k,2)-node(midnode)%value_c(k,1))*t1_tip(k)
                        delta_u_I_endnode =delta_u_I_endnode +(node(endnode)%value_c(k,2)-node(endnode)%value_c(k,1))*n_tip(k)
                        delta_u_II_endnode=delta_u_II_endnode+(node(endnode)%value_c(k,2)-node(endnode)%value_c(k,1))*t1_tip(k)
                      end do
                      ! Calculate K_I and K_II using extrapolation of K_I and K_II at mid node and end node.
                      !K_I =mu/(kappa+1.0d0)*dsqrt(c_2pi)/(ratio-1.d0)*(ratio*(delta_u_I_endnode )/dsqrt(r_endnode)-(delta_u_I_midnode )/dsqrt(r_midnode))
                      !K_II=mu/(kappa+1.0d0)*dsqrt(c_2pi)/(ratio-1.d0)*(ratio*(delta_u_II_endnode)/dsqrt(r_endnode)-(delta_u_II_midnode)/dsqrt(r_midnode))
                      ! Calculate K_I and K_II using the end node solution
                      !K_I =mu/(kappa+1.0d0)*dsqrt(c_2pi/r_endnode)*delta_u_I_endnode
                      !K_II=mu/(kappa+1.0d0)*dsqrt(c_2pi/r_endnode)*delta_u_II_endnode
                      ! Calculate K_I and K_II using the mid node solution
                      K_I =mu/(kappa+1.0d0)*sqrt(c_2pi)*delta_u_I_midnode/dsqrt(r_midnode)
                      K_II=mu/(kappa+1.0d0)*sqrt(c_2pi)*delta_u_II_midnode/dsqrt(r_midnode)
                      ! Using the comparison between analytical and numerical (modified quarter-point)
                      !K_I =mu/(kappa+1.0d0)*sqrt(c_2pi/r_endnode)*(-1.d0/3.d0*delta_u_I_endnode +8.d0/3.d0*delta_u_I_midnode )
                      !K_II=mu/(kappa+1.0d0)*sqrt(c_2pi/r_endnode)*(-1.d0/3.d0*delta_u_II_endnode+8.d0/3.d0*delta_u_II_midnode)
                      ! Write results
                      write(output_fileunit,fmtstr) node(kn)%id, omega, zabs(K_I), zabs(K_II), fbem_zarg(K_I), fbem_zarg(K_II)

                      !
                      ! Esto es temporal!!!
                      !
                      !write(output_fileunit,'(2e25.16)') omega, zabs(delta_u_I_midnode)/dsqrt(r_midnode)

                    !
                    ! Ordinary boundary
                    !
                    case (fbem_boundary_class_ordinary)

                      ! Calculate t tangential and normal
                      t_I=0.d0
                      t_II=0.d0
                      do k=1,problem%n
                        t_I =t_I +node(kn)%value_c(problem%n+1+k,1)*n_tip(k)
                        t_II=t_II+node(kn)%value_c(problem%n+1+k,1)*t1_tip(k)
                      end do
                      ! Calculate K_I and K_II using the mid node solution
                      K_I =sqrt(c_2pi*r_endnode)*(t_I)
                      K_II=sqrt(c_2pi*r_endnode)*(t_II)
                      ! Write results
                      write(output_fileunit,fmtstr) node(kn)%id, omega, zabs(K_I), zabs(K_II), fbem_zarg(K_I), fbem_zarg(K_II)

                  end select

          end select

        !
        ! Other elements
        !
        case default
          call fbem_error_message(error_unit,0,'element',element(se)%id,&
                                  'this type of element of the node does not admit a SIF calculation')
      end select
    end if
  end do

end subroutine export_solution_mechanics_harmonic_sif
