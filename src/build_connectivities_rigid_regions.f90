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

subroutine build_connectivities_rigid_regions

  ! Fortran 2003 intrinsic module
  use iso_fortran_env

  ! fbem modules
  use fbem_numerical
  use fbem_geometry
  use fbem_string_handling
  use fbem_shape_functions
  use fbem_symmetry
  use fbem_data_structures

  ! Problem variables module
  use problem_variables

  ! No implicit variables are allowed
  implicit none

  ! Local variables
  integer                               :: sn, sp
  integer                               :: kr, kn, kp

  if (verbose_level.ge.2) call fbem_timestamp_w_message(output_unit,2,'START building connectivities in rigid regions')

  do kr=1,n_regions
    if ((region(kr)%class.eq.fbem_fe).and.(region(kr)%type.eq.fbem_rigid)) then
      ! Check if the any node of the region is a master node already, which is not allowed.
      do kp=1,region(kr)%n_fe_subregions
        sp=fe_subregion(region(kr)%fe_subregion(kp))%part
        do kn=1,part(sp)%n_nodes
          sn=part(sp)%node(kn)
          if (node(sn)%rigid_link.ne.0) then
            call fbem_error_message(error_unit,0,'region',region(kr)%id,'direct connection between rigid regions is not allowed')
          end if
        end do
      end do
      ! If the master node id is known, copy its coordinates
      if (region(kr)%master_node.ne.0) then
        allocate (region(kr)%master_node_x(problem%n))
        region(kr)%master_node_x=node(region(kr)%master_node)%x
      end if
      ! Find the master node, which must not share position with any other node.
      do kp=1,region(kr)%n_fe_subregions
        sp=fe_subregion(region(kr)%fe_subregion(kp))%part
        do kn=1,part(sp)%n_nodes
          sn=part(sp)%node(kn)
          if (sqrt(sum((node(sn)%x-region(kr)%master_node_x)**2)).lt.geometric_tolerance) then
            if ((region(kr)%master_node.eq.0).or.(region(kr)%master_node.eq.sn)) then
              node(sn)%rigid_link=1
              region(kr)%master_node=sn
            else
              call fbem_error_message(error_unit,0,'region',region(kr)%id,'at the master node position there are more than one node')
            end if
          end if
        end do
      end do

      !write(*,*) 'NODO MAESTRO'
      !write(*,*) region(kr)%master_node
      !write(*,*) region(kr)%master_node_x

      ! Copy the master node iid to all subregions
      do kp=1,region(kr)%n_fe_subregions
        fe_subregion(region(kr)%fe_subregion(kp))%master_node=region(kr)%master_node
      end do
      ! Count the number of slave nodes,
      node(region(kr)%master_node)%n_slaves=0
      do kp=1,region(kr)%n_fe_subregions
        sp=fe_subregion(region(kr)%fe_subregion(kp))%part
        do kn=1,part(sp)%n_nodes
          sn=part(sp)%node(kn)
          if (node(sn)%rigid_link.eq.0) then
            node(region(kr)%master_node)%n_slaves=node(region(kr)%master_node)%n_slaves+1
            node(sn)%rigid_link=2
          end if
        end do
      end do
      ! Allocate slave nodes vector in master node
      allocate(node(region(kr)%master_node)%slave(node(region(kr)%master_node)%n_slaves))
      ! Reset rigid_link in slave nodes
      do kp=1,region(kr)%n_fe_subregions
        sp=fe_subregion(region(kr)%fe_subregion(kp))%part
        do kn=1,part(sp)%n_nodes
          sn=part(sp)%node(kn)
          if (node(sn)%rigid_link.eq.2) then
            node(sn)%rigid_link=0
          end if
        end do
      end do
      ! Build the slave nodes vector in master node and the master node in each slave node.
      node(region(kr)%master_node)%n_slaves=0
      do kp=1,region(kr)%n_fe_subregions
        sp=fe_subregion(region(kr)%fe_subregion(kp))%part
        do kn=1,part(sp)%n_nodes
          sn=part(sp)%node(kn)
          if (node(sn)%rigid_link.eq.0) then
            node(region(kr)%master_node)%n_slaves=node(region(kr)%master_node)%n_slaves+1
            node(region(kr)%master_node)%slave(node(region(kr)%master_node)%n_slaves)=sn
            node(sn)%rigid_link=2
            node(sn)%master=region(kr)%master_node
            node(sn)%n_dof=3*(problem%n-1) ! all slave nodes, if they have rotations, the must be global rotations (especially shells)
          end if
        end do
      end do
    end if
  end do

  if (verbose_level.ge.2) call fbem_timestamp_w_message(output_unit,2,'END building connectivities in rigid regions')


end subroutine build_connectivities_rigid_regions
