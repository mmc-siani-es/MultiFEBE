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

!! @author Jacob David Rodriguez Bordon (jacobdavid.rodriguezbordon@ulpgc.es)
!!
!! @version 2.0
!!
!! <b>This module implements beam, bar, spring/dashpot finite elements.</b>
module fbem_fem_beams

  ! Fortran 2003 standard
  use iso_fortran_env

  ! fbem modules
  use fbem_numerical
  use fbem_string_handling
  use fbem_quad_rules
  use fbem_shape_functions
  use fbem_geometry

  implicit none

  private

  ! ================================================================================================================================
  ! DISCRETE TRANSLATIONAL SPRING ELEMENT (DISTRA)
  !
  ! TYPE OF MESH ELEMENT:
  !
  ! DOF per node:
  !   2D: (u1,u2)'
  !   3D: (u1,u2,u3)'
  !
  public :: fbem_fem_distra_L
  public :: fbem_fem_distra_K_real
  public :: fbem_fem_distra_K_complex
  public :: fbem_fem_distra_K_static
  public :: fbem_fem_distra_K_harmonic
  public :: fbem_fem_distra_KLT_static
  public :: fbem_fem_distra_KLT_harmonic
  !
  !public :: fbem_fem_distra_test
  ! ================================================================================================================================

  ! ================================================================================================================================
  ! DISCRETE ROTATIONAL/TRANSLATIONAL SPRING ELEMENT (DISROTRA)
  !
  ! TYPE OF MESH ELEMENT:
  !
  ! LINE2: (1)---------------(2)
  !
  ! DOF per node:
  !   2D: (u1,u2,theta3)'
  !   3D: (u1,u2,u3,theta1,theta2,theta3)'
  !
  public :: fbem_fem_disrotra_L
  public :: fbem_fem_disrotra_K_real
  public :: fbem_fem_disrotra_K_complex
  public :: fbem_fem_disrotra_K_static
  public :: fbem_fem_disrotra_K_harmonic
  public :: fbem_fem_disrotra_KLT_static
  public :: fbem_fem_disrotra_KLT_harmonic
  !
  !public :: fbem_fem_disrotra_test
  ! ================================================================================================================================

  ! ================================================================================================================================
  ! STRAIGHT BAR ELEMENT (BAR)
  !
  ! TYPE OF MESH ELEMENT:
  !
  ! LINE2: (1)---------------(2)
  !
  ! DOF per node:
  !   2D: (u1,u2)'
  !   3D: (u1,u2,u3)'
  !
  public :: fbem_fem_bar_L
  public :: fbem_fem_bar_K
  public :: fbem_fem_bar_M
  public :: fbem_fem_bar_K_static
  public :: fbem_fem_bar_K_harmonic
  public :: fbem_fem_bar_DB
  public :: fbem_fem_bar_stress_resultants
  !
  !public :: fbem_fem_bar_test
  ! ================================================================================================================================

  ! ================================================================================================================================
  ! STRAIGHT EULER-BERNOULLI AND TIMOSHENKO BEAM ELEMENTS (STRBEAM)
  !
  ! Type of elements (straight):
  !
  ! LINE2: (1)---------------(2)
  !
  ! LINE3: (1)------(3)------(2)
  !
  ! DOF per node:
  !   2D: (u1,u2,theta3)'
  !   3D: (u1,u2,u3,theta1,theta2,theta3)'
  ! DOF for mid-node of line3 with subtype=0:
  !   2D: (u1,u2)'
  !   3D: (u1,u2,u3)'
  ! DOF for mid-nodes of line4 with subtype=0:
  !   2D: (u1')'
  !
  public :: fbem_fem_strbeam_L_element
  public :: fbem_fem_strbeam_Na
  public :: fbem_fem_strbeam_Ka
  public :: fbem_fem_strbeam_Ma
  public :: fbem_fem_strbeam_Qa
  public :: fbem_fem_strbeam_Nl
  public :: fbem_fem_strbeam_Kl
  public :: fbem_fem_strbeam_Ml
  public :: fbem_fem_strbeam_Ql
  public :: fbem_fem_strbeam_K
  public :: fbem_fem_strbeam_M
  public :: fbem_fem_strbeam_Madd
  public :: fbem_fem_strbeam_Q
  public :: fbem_fem_strbeam_DBa
  public :: fbem_fem_strbeam_DBl
  public :: fbem_fem_strbeam_DB
  !
  public :: fbem_fem_strbeam_N
  public :: fbem_fem_strbeam_K_static
  public :: fbem_fem_strbeam_K_harmonic
  public :: fbem_fem_strbeam_Q_midline
  public :: fbem_fem_strbeam_L_load
  public :: fbem_fem_strbeam_stress_resultants
  !
  !public :: fbem_fem_strbeam_test
  ! ================================================================================================================================

  ! ================================================================================================================================
  ! DEGENERATED BEAM ELEMENTS (DEGBEAM)
  !
  ! Type of elements (straight or curved):
  !
  ! LINE2: (1)---------------(2)
  !
  ! LINE3: (1)------(3)------(2)
  !
  ! LINE4: (1)---(3)---(4)---(2)
  !
  ! DOF per node:
  !   2D: (u1,u2,theta3)'
  !   3D: (u1,u2,u3,theta1,theta2,theta3)'
  !
  public :: fbem_fem_degbeam3d_x
  public :: fbem_fem_degbeam3d_L_element
  public :: fbem_fem_degbeam3d_L_load
  public :: fbem_fem_degbeam3d_K
  public :: fbem_fem_degbeam3d_Q_body
  public :: fbem_fem_degbeam3d_Q_midline
  public :: fbem_fem_degbeam3d_stress_resultants
  public :: fbem_fem_degbeam3d_stress_resultants_ep
  !
  public :: fbem_fem_degbeam_K_static
  public :: fbem_fem_degbeam_K_harmonic
  public :: fbem_fem_degbeam_Q_midline
  public :: fbem_fem_degbeam_Q_body
  public :: fbem_fem_degbeam_L_load
  !
  public :: fbem_fem_degbeam_stress_resultants
  ! ================================================================================================================================

contains

  ! ================================================================================================================================
  ! DISCRETE TRANSLATIONAL SPRING ELEMENT (DISTRA): 2-node structural element with only displacements DOFs.
  ! ================================================================================================================================

  !! Coordinate transformation matrix L
  function fbem_fem_distra_L(rn,local_axes,nodal_axes)
    implicit none
    integer           :: rn                                    !! Number of dimensions of the ambient space
    real(kind=real64) :: local_axes(rn,rn)                     !! Local axes where springs and dashpots are defined
    real(kind=real64) :: nodal_axes(rn,rn,2)                   !! Nodal axes for each node DOFs nodal_axes(component,axis,node)
    real(kind=real64) :: fbem_fem_distra_L(2*rn,2*rn)
    integer           :: kn
    real(kind=real64) :: e(rn,rn), Ln(rn,rn)
    fbem_fem_distra_L=0.d0
    do kn=1,2
      e=nodal_axes(:,:,kn)
      call fbem_coordinate_transformation_L(rn,local_axes,e,Ln)
      fbem_fem_distra_L((rn*(kn-1)+1):(rn*kn),(rn*(kn-1)+1):(rn*kn))=Ln
    end do
  end function fbem_fem_distra_L

  !! Stiffness matrix in local axes
  function fbem_fem_distra_K_real(rn,K)
    implicit none
    integer           :: rn
    real(kind=real64) :: K(8)                              !! Springs stiffnesses
    real(kind=real64) :: fbem_fem_distra_K_real(2*rn,2*rn) !! Stiffness matrix in local coordinates
    fbem_fem_distra_K_real=0.d0
    select case (rn)
      case (2)
        ! Axis 1'
        fbem_fem_distra_K_real(1,1)= K(1)
        fbem_fem_distra_K_real(1,3)=-K(1)
        fbem_fem_distra_K_real(3,1)=-K(1)
        fbem_fem_distra_K_real(3,3)= K(1)
        ! Axis 2'
        fbem_fem_distra_K_real(2,2)= K(2)
        fbem_fem_distra_K_real(2,4)=-K(2)
        fbem_fem_distra_K_real(4,2)=-K(2)
        fbem_fem_distra_K_real(4,4)= K(2)
      case (3)
        ! Axis 1'
        fbem_fem_distra_K_real(1,1)= K(1)
        fbem_fem_distra_K_real(1,4)=-K(1)
        fbem_fem_distra_K_real(4,1)=-K(1)
        fbem_fem_distra_K_real(4,4)= K(1)
        ! Axis 2'
        fbem_fem_distra_K_real(2,2)= K(2)
        fbem_fem_distra_K_real(2,5)=-K(2)
        fbem_fem_distra_K_real(5,2)=-K(2)
        fbem_fem_distra_K_real(5,5)= K(2)
        ! Axis 3'
        fbem_fem_distra_K_real(3,3)= K(3)
        fbem_fem_distra_K_real(3,6)=-K(3)
        fbem_fem_distra_K_real(6,3)=-K(3)
        fbem_fem_distra_K_real(6,6)= K(3)
    end select
  end function fbem_fem_distra_K_real

  !! Stiffness matrix in local axes
  function fbem_fem_distra_K_complex(rn,K)
    implicit none
    integer              :: rn
    complex(kind=real64) :: K(8)                                 !! Springs copmlex stiffnesses
    complex(kind=real64) :: fbem_fem_distra_K_complex(2*rn,2*rn) !! Stiffness matrix in local coordinates
    fbem_fem_distra_K_complex=0.d0
    select case (rn)
      case (2)
        ! Axis 1'
        fbem_fem_distra_K_complex(1,1)= K(1)
        fbem_fem_distra_K_complex(1,3)=-K(1)
        fbem_fem_distra_K_complex(3,1)=-K(1)
        fbem_fem_distra_K_complex(3,3)= K(1)
        ! Axis 2'
        fbem_fem_distra_K_complex(2,2)= K(2)
        fbem_fem_distra_K_complex(2,4)=-K(2)
        fbem_fem_distra_K_complex(4,2)=-K(2)
        fbem_fem_distra_K_complex(4,4)= K(2)
      case (3)
        ! Axis 1'
        fbem_fem_distra_K_complex(1,1)= K(1)
        fbem_fem_distra_K_complex(1,4)=-K(1)
        fbem_fem_distra_K_complex(4,1)=-K(1)
        fbem_fem_distra_K_complex(4,4)= K(1)
        ! Axis 2'
        fbem_fem_distra_K_complex(2,2)= K(2)
        fbem_fem_distra_K_complex(2,5)=-K(2)
        fbem_fem_distra_K_complex(5,2)=-K(2)
        fbem_fem_distra_K_complex(5,5)= K(2)
        ! Axis 3'
        fbem_fem_distra_K_complex(3,3)= K(3)
        fbem_fem_distra_K_complex(3,6)=-K(3)
        fbem_fem_distra_K_complex(6,3)=-K(3)
        fbem_fem_distra_K_complex(6,6)= K(3)
    end select
  end function fbem_fem_distra_K_complex

  !! Calculation of all matrices needed for static analysis
  subroutine fbem_fem_distra_K_static(rn,local_axis,nodal_axes,Kval,K)
    implicit none
    ! I/O
    integer           :: rn                  !! Number of dimensions of the ambient space
    real(kind=real64) :: local_axis(rn,rn)   !! Local axes where springs and dashpots are defined
    real(kind=real64) :: nodal_axes(rn,rn,2) !! Axes for each node DOFs nodal_axes(component,axis,node)
    real(kind=real64) :: Kval(8)             !! Springs stiffnesses
    real(kind=real64) :: K(2*rn,2*rn)        !! Stiffness matrix
    ! Local
    real(kind=real64) :: L(2*rn,2*rn)
    real(kind=real64) :: Kl(2*rn,2*rn)
    ! Calculation of matrices
    L=fbem_fem_distra_L(rn,local_axis,nodal_axes)
    Kl=fbem_fem_distra_K_real(rn,Kval)
    K=matmul(L,matmul(Kl,transpose(L)))
  end subroutine fbem_fem_distra_K_static

  !! Calculation of all matrices needed for time harmonic analysis
  subroutine fbem_fem_distra_K_harmonic(rn,local_axis,nodal_axes,Kval,K)
    implicit none
    ! I/O
    integer              :: rn                  !! Number of dimensions of the ambient space
    real(kind=real64)    :: local_axis(rn,rn)   !! Local axes where springs and dashpots are defined
    real(kind=real64)    :: nodal_axes(rn,rn,2) !! Axes for each node DOFs nodal_axes(component,axis,node)
    complex(kind=real64) :: Kval(8)             !! Springs complex stiffnesses
    complex(kind=real64) :: K(2*rn,2*rn)        !! Complex stiffness matrix
    ! Local
    real(kind=real64)    :: L(2*rn,2*rn)
    complex(kind=real64) :: Kl(2*rn,2*rn)
    ! Calculation of matrices
    L=fbem_fem_distra_L(rn,local_axis,nodal_axes)
    Kl=fbem_fem_distra_K_complex(rn,Kval)
    K=matmul(L,matmul(Kl,transpose(L)))
  end subroutine fbem_fem_distra_K_harmonic

  !! Calculation of matrix for spring force calculation
  subroutine fbem_fem_distra_KLT_static(rn,local_axis,nodal_axes,Kval,KLT)
    implicit none
    ! I/O
    integer           :: rn                  !! Number of dimensions of the ambient space
    real(kind=real64) :: local_axis(rn,rn)   !! Local axes where springs and dashpots are defined
    real(kind=real64) :: nodal_axes(rn,rn,2) !! Axes for each node DOFs nodal_axes(component,axis,node)
    real(kind=real64) :: Kval(8)             !! Springs stiffnesses
    real(kind=real64) :: KLT(2*rn,2*rn)        !! Stiffness matrix
    ! Local
    real(kind=real64) :: L(2*rn,2*rn)
    real(kind=real64) :: Kl(2*rn,2*rn)
    ! Calculation of matrices
    L=fbem_fem_distra_L(rn,local_axis,nodal_axes)
    Kl=fbem_fem_distra_K_real(rn,Kval)
    KLT=matmul(Kl,transpose(L))
  end subroutine fbem_fem_distra_KLT_static

  !! Calculation of matrix for spring force calculation
  subroutine fbem_fem_distra_KLT_harmonic(rn,local_axis,nodal_axes,Kval,KLT)
    implicit none
    ! I/O
    integer              :: rn                  !! Number of dimensions of the ambient space
    real(kind=real64)    :: local_axis(rn,rn)   !! Local axes where springs and dashpots are defined
    real(kind=real64)    :: nodal_axes(rn,rn,2) !! Axes for each node DOFs nodal_axes(component,axis,node)
    complex(kind=real64) :: Kval(8)             !! Springs stiffnesses
    complex(kind=real64) :: KLT(2*rn,2*rn)      !! Stiffness matrix
    ! Local
    real(kind=real64) :: L(2*rn,2*rn)
    complex(kind=real64) :: Kl(2*rn,2*rn)
    ! Calculation of matrices
    L=fbem_fem_distra_L(rn,local_axis,nodal_axes)
    Kl=fbem_fem_distra_K_complex(rn,Kval)
    KLT=matmul(Kl,transpose(L))
  end subroutine fbem_fem_distra_KLT_harmonic

!  ! Tests
!  subroutine fbem_fem_distra_test()
!    implicit none
!    integer, parameter   :: rn=3
!    real(kind=real64)    :: local_axes(rn,rn)   !! Local axes where springs and dashpots are defined
!    real(kind=real64)    :: nodal_axes(rn,rn,2) !! Axes for each node DOFs nodal_axes(component,axis,node)
!    real(kind=real64)    :: Kval(rn)
!    real(kind=real64)    :: K(2*rn,2*rn)        !! Complex stiffness matrix
!    integer              :: kn
!    real(kind=real64)    :: angle
!    ! Static test
!    angle=c_pi/2.d0
!    local_axes(:,1)=[cos(angle),sin(angle),0.d0]
!    local_axes(:,2)=[-sin(angle),cos(angle),0.d0]
!    local_axes(:,3)=[0.d0,0.d0,1.d0]
!    nodal_axes(:,1,1)=local_axes(:,1)
!    nodal_axes(:,2,1)=local_axes(:,2)
!    nodal_axes(:,3,1)=local_axes(:,3)
!    nodal_axes(:,1,2)=local_axes(:,1)
!    nodal_axes(:,2,2)=local_axes(:,2)
!    nodal_axes(:,3,2)=local_axes(:,3)
!!    nodal_axes(:,1,1)=[1,0,0]
!!    nodal_axes(:,2,1)=[0,1,0]
!!    nodal_axes(:,3,1)=[0,0,1]
!!    nodal_axes(:,1,2)=[1,0,0]
!!    nodal_axes(:,2,2)=[0,1,0]
!!    nodal_axes(:,3,2)=[0,0,1]
!    Kval=[3,2,1]
!    call fbem_fem_distra_K_static(rn,local_axes,nodal_axes,Kval,K)
!    write(*,*) 'K'
!    do kn=1,2*rn
!      write(*,*) K(kn,:)
!    end do
!  end subroutine fbem_fem_distra_test

  ! ================================================================================================================================
  ! DISCRETE ROTATIONAL/TRANSLATIONAL SPRING ELEMENT (DISROTRA): 2-node structural element with displacements and rotational DOFs.
  ! ================================================================================================================================

  !! Coordinate transformation matrix L
  function fbem_fem_disrotra_L(rn,local_axes,nodal_axes)
    implicit none
    integer           :: rn                                         !! Number of dimensions of the ambient space
    real(kind=real64) :: local_axes(rn,rn)                          !! Local axes where springs and dashpots are defined
    real(kind=real64) :: nodal_axes(rn,rn,2)                        !! Nodal axes for each node DOFs nodal_axes(component,axis,node)
    real(kind=real64) :: fbem_fem_disrotra_L(2*3*(rn-1),2*3*(rn-1))
    real(kind=real64) :: e(rn,rn), Ln(rn,rn)
    fbem_fem_disrotra_L=0.d0
    select case (rn)
    case (2)
      e=nodal_axes(:,:,1)
      call fbem_coordinate_transformation_L(2,local_axes,e,Ln)
      fbem_fem_disrotra_L(1:2,1:2)=Ln
      fbem_fem_disrotra_L(  3,  3)=1.d0
      e=nodal_axes(:,:,2)
      call fbem_coordinate_transformation_L(2,local_axes,e,Ln)
      fbem_fem_disrotra_L(4:5,4:5)=Ln
      fbem_fem_disrotra_L(  6,  6)=1.d0
    case (3)
      e=nodal_axes(:,:,1)
      call fbem_coordinate_transformation_L(3,local_axes,e,Ln)
      fbem_fem_disrotra_L(1:3,1:3)=Ln
      fbem_fem_disrotra_L(4:6,4:6)=Ln
      e=nodal_axes(:,:,2)
      call fbem_coordinate_transformation_L(3,local_axes,e,Ln)
      fbem_fem_disrotra_L(7:9,7:9)=Ln
      fbem_fem_disrotra_L(10:12,10:12)=Ln
    end select
  end function fbem_fem_disrotra_L

  !! Stiffness matrix in local axes
  function fbem_fem_disrotra_K_real(rn,K)
    implicit none
    integer           :: rn
    real(kind=real64) :: K(8)                                            !! Springs stiffnesses
    real(kind=real64) :: fbem_fem_disrotra_K_real(2*3*(rn-1),2*3*(rn-1)) !! Stiffness matrix in local coordinates
    fbem_fem_disrotra_K_real=0.d0
    select case (rn)
      case (2)
        ! Translation axis 1'
        fbem_fem_disrotra_K_real(1,1)= K(1)
        fbem_fem_disrotra_K_real(1,4)=-K(1)
        fbem_fem_disrotra_K_real(4,1)=-K(1)
        fbem_fem_disrotra_K_real(4,4)= K(1)
        ! Coupled translation axis 2' / rotation axis 3
        fbem_fem_disrotra_K_real(2,2)= K(2)
        fbem_fem_disrotra_K_real(2,3)= K(4)
        fbem_fem_disrotra_K_real(3,2)= K(4)
        fbem_fem_disrotra_K_real(2,5)=-K(2)
        fbem_fem_disrotra_K_real(5,2)=-K(2)
        fbem_fem_disrotra_K_real(2,6)= K(4)
        fbem_fem_disrotra_K_real(6,2)= K(4)
        fbem_fem_disrotra_K_real(3,3)= K(3)
        fbem_fem_disrotra_K_real(3,5)=-K(4)
        fbem_fem_disrotra_K_real(5,3)=-K(4)
        fbem_fem_disrotra_K_real(3,6)=2.d0*K(4)**2/K(2)-K(3)
        fbem_fem_disrotra_K_real(6,3)=2.d0*K(4)**2/K(2)-K(3)
        fbem_fem_disrotra_K_real(5,5)= K(2)
        fbem_fem_disrotra_K_real(5,6)=-K(4)
        fbem_fem_disrotra_K_real(6,5)=-K(4)
        fbem_fem_disrotra_K_real(6,6)= K(3)
      case (3)
        ! Translation axis 1'
        fbem_fem_disrotra_K_real( 1, 1)= K(1)
        fbem_fem_disrotra_K_real( 1, 7)=-K(1)
        fbem_fem_disrotra_K_real( 7, 1)=-K(1)
        fbem_fem_disrotra_K_real( 7, 7)= K(1)
        ! Rotation axis 1'
        fbem_fem_disrotra_K_real( 4, 4)= K(4)
        fbem_fem_disrotra_K_real( 4,10)=-K(4)
        fbem_fem_disrotra_K_real(10, 4)=-K(4)
        fbem_fem_disrotra_K_real(10,10)= K(4)
        ! Coupled translation axis 2' / rotation axis 3'
        fbem_fem_disrotra_K_real( 2, 2)= K(2)
        fbem_fem_disrotra_K_real( 2, 6)= K(7)
        fbem_fem_disrotra_K_real( 6, 2)= K(7)
        fbem_fem_disrotra_K_real( 2, 8)=-K(2)
        fbem_fem_disrotra_K_real( 8, 2)=-K(2)
        fbem_fem_disrotra_K_real( 2,12)= K(7)
        fbem_fem_disrotra_K_real(12, 2)= K(7)
        fbem_fem_disrotra_K_real( 6, 6)= K(6)
        fbem_fem_disrotra_K_real( 6, 8)=-K(7)
        fbem_fem_disrotra_K_real( 8, 6)=-K(7)
        fbem_fem_disrotra_K_real( 6,12)=2.d0*K(7)**2/K(2)-K(6)
        fbem_fem_disrotra_K_real(12, 6)=2.d0*K(7)**2/K(2)-K(6)
        fbem_fem_disrotra_K_real( 8, 8)= K(2)
        fbem_fem_disrotra_K_real( 8,12)=-K(7)
        fbem_fem_disrotra_K_real(12, 8)=-K(7)
        fbem_fem_disrotra_K_real(12,12)= K(6)
        ! Coupled translation axis 3' / rotation axis 2'
        fbem_fem_disrotra_K_real( 3, 3)= K(3)
        fbem_fem_disrotra_K_real( 3, 5)=-K(8)
        fbem_fem_disrotra_K_real( 5, 3)=-K(8)
        fbem_fem_disrotra_K_real( 3, 9)=-K(3)
        fbem_fem_disrotra_K_real( 9, 3)=-K(3)
        fbem_fem_disrotra_K_real( 3,11)=-K(8)
        fbem_fem_disrotra_K_real(11, 3)=-K(8)
        fbem_fem_disrotra_K_real( 5, 5)= K(5)
        fbem_fem_disrotra_K_real( 5, 9)= K(8)
        fbem_fem_disrotra_K_real( 9, 5)= K(8)
        fbem_fem_disrotra_K_real( 5,11)=2.d0*K(8)**2/K(3)-K(5)
        fbem_fem_disrotra_K_real(11, 5)=2.d0*K(8)**2/K(3)-K(5)
        fbem_fem_disrotra_K_real( 9, 9)= K(3)
        fbem_fem_disrotra_K_real( 9,11)= K(8)
        fbem_fem_disrotra_K_real(11, 9)= K(8)
        fbem_fem_disrotra_K_real(11,11)= K(5)
    end select
  end function fbem_fem_disrotra_K_real

  !! Stiffness matrix in local axes
  function fbem_fem_disrotra_K_complex(rn,K)
    implicit none
    integer              :: rn
    complex(kind=real64) :: K(8)                                               !! Springs stiffnesses
    complex(kind=real64) :: fbem_fem_disrotra_K_complex(2*3*(rn-1),2*3*(rn-1)) !! Stiffness matrix in local coordinates
    fbem_fem_disrotra_K_complex=0.d0
    select case (rn)
      case (2)
        ! Translation axis 1'
        fbem_fem_disrotra_K_complex(1,1)= K(1)
        fbem_fem_disrotra_K_complex(1,4)=-K(1)
        fbem_fem_disrotra_K_complex(4,1)=-K(1)
        fbem_fem_disrotra_K_complex(4,4)= K(1)
        ! Coupled translation axis 2' / rotation axis 3
        fbem_fem_disrotra_K_complex(2,2)= K(2)
        fbem_fem_disrotra_K_complex(2,3)= K(4)
        fbem_fem_disrotra_K_complex(3,2)= K(4)
        fbem_fem_disrotra_K_complex(2,5)=-K(2)
        fbem_fem_disrotra_K_complex(5,2)=-K(2)
        fbem_fem_disrotra_K_complex(2,6)= K(4)
        fbem_fem_disrotra_K_complex(6,2)= K(4)
        fbem_fem_disrotra_K_complex(3,3)= K(3)
        fbem_fem_disrotra_K_complex(3,5)=-K(4)
        fbem_fem_disrotra_K_complex(5,3)=-K(4)
        fbem_fem_disrotra_K_complex(3,6)=2.d0*K(4)**2/K(2)-K(3)
        fbem_fem_disrotra_K_complex(6,3)=2.d0*K(4)**2/K(2)-K(3)
        fbem_fem_disrotra_K_complex(5,5)= K(2)
        fbem_fem_disrotra_K_complex(5,6)=-K(4)
        fbem_fem_disrotra_K_complex(6,5)=-K(4)
        fbem_fem_disrotra_K_complex(6,6)= K(3)
      case (3)
        ! Translation axis 1'
        fbem_fem_disrotra_K_complex( 1, 1)= K(1)
        fbem_fem_disrotra_K_complex( 1, 7)=-K(1)
        fbem_fem_disrotra_K_complex( 7, 1)=-K(1)
        fbem_fem_disrotra_K_complex( 7, 7)= K(1)
        ! Rotation axis 1'
        fbem_fem_disrotra_K_complex( 4, 4)= K(4)
        fbem_fem_disrotra_K_complex( 4,10)=-K(4)
        fbem_fem_disrotra_K_complex(10, 4)=-K(4)
        fbem_fem_disrotra_K_complex(10,10)= K(4)
        ! Coupled translation axis 2' / rotation axis 3'
        fbem_fem_disrotra_K_complex( 2, 2)= K(2)
        fbem_fem_disrotra_K_complex( 2, 6)= K(7)
        fbem_fem_disrotra_K_complex( 6, 2)= K(7)
        fbem_fem_disrotra_K_complex( 2, 8)=-K(2)
        fbem_fem_disrotra_K_complex( 8, 2)=-K(2)
        fbem_fem_disrotra_K_complex( 2,12)= K(7)
        fbem_fem_disrotra_K_complex(12, 2)= K(7)
        fbem_fem_disrotra_K_complex( 6, 6)= K(6)
        fbem_fem_disrotra_K_complex( 6, 8)=-K(7)
        fbem_fem_disrotra_K_complex( 8, 6)=-K(7)
        fbem_fem_disrotra_K_complex( 6,12)=2.d0*K(7)**2/K(2)-K(6)
        fbem_fem_disrotra_K_complex(12, 6)=2.d0*K(7)**2/K(2)-K(6)
        fbem_fem_disrotra_K_complex( 8, 8)= K(2)
        fbem_fem_disrotra_K_complex( 8,12)=-K(7)
        fbem_fem_disrotra_K_complex(12, 8)=-K(7)
        fbem_fem_disrotra_K_complex(12,12)= K(6)
        ! Coupled translation axis 3' / rotation axis 2'
        fbem_fem_disrotra_K_complex( 3, 3)= K(3)
        fbem_fem_disrotra_K_complex( 3, 5)=-K(8)
        fbem_fem_disrotra_K_complex( 5, 3)=-K(8)
        fbem_fem_disrotra_K_complex( 3, 9)=-K(3)
        fbem_fem_disrotra_K_complex( 9, 3)=-K(3)
        fbem_fem_disrotra_K_complex( 3,11)=-K(8)
        fbem_fem_disrotra_K_complex(11, 3)=-K(8)
        fbem_fem_disrotra_K_complex( 5, 5)= K(5)
        fbem_fem_disrotra_K_complex( 5, 9)= K(8)
        fbem_fem_disrotra_K_complex( 9, 5)= K(8)
        fbem_fem_disrotra_K_complex( 5,11)=2.d0*K(8)**2/K(3)-K(5)
        fbem_fem_disrotra_K_complex(11, 5)=2.d0*K(8)**2/K(3)-K(5)
        fbem_fem_disrotra_K_complex( 9, 9)= K(3)
        fbem_fem_disrotra_K_complex( 9,11)= K(8)
        fbem_fem_disrotra_K_complex(11, 9)= K(8)
        fbem_fem_disrotra_K_complex(11,11)= K(5)
    end select
  end function fbem_fem_disrotra_K_complex

  !! Calculation of all matrices needed for static analysis
  subroutine fbem_fem_disrotra_K_static(rn,local_axis,nodal_axes,Kval,K)
    implicit none
    ! I/O
    integer           :: rn                       !! Number of dimensions of the ambient space
    real(kind=real64) :: local_axis(rn,rn)        !! Local axes where springs and dashpots are defined
    real(kind=real64) :: nodal_axes(rn,rn,2)      !! Axes for each node DOFs nodal_axes(component,axis,node)
    real(kind=real64) :: Kval(8)                  !! Springs stiffnesses: 2D (k1,k2,kr3,k2r3), 3D (k1,k2,k3,kr1,kr2,kr3,k2r3,k3r2)
    real(kind=real64) :: K(2*3*(rn-1),2*3*(rn-1)) !! Stiffness matrix
    ! Local
    real(kind=real64) :: L(2*3*(rn-1),2*3*(rn-1))
    real(kind=real64) :: Kl(2*3*(rn-1),2*3*(rn-1))
    ! Calculation of matrices
    L=fbem_fem_disrotra_L(rn,local_axis,nodal_axes)
    Kl=fbem_fem_disrotra_K_real(rn,Kval)
    K=matmul(L,matmul(Kl,transpose(L)))
  end subroutine fbem_fem_disrotra_K_static

  !! Calculation of all matrices needed for time harmonic analysis
  subroutine fbem_fem_disrotra_K_harmonic(rn,local_axis,nodal_axes,Kval,K)
    implicit none
    ! I/O
    integer              :: rn                       !! Number of dimensions of the ambient space
    real(kind=real64)    :: local_axis(rn,rn)        !! Local axes where springs and dashpots are defined
    real(kind=real64)    :: nodal_axes(rn,rn,2)      !! Axes for each node DOFs nodal_axes(component,axis,node)
    complex(kind=real64) :: Kval(8)                  !! Springs stiffnesses: 2D (k1,k2,kr3,k2r3), 3D (k1,k2,k3,kr1,kr2,kr3,k2r3,k3r2)
    complex(kind=real64) :: K(2*3*(rn-1),2*3*(rn-1)) !! Complex stiffness matrix
    ! Local
    real(kind=real64)    :: L(2*3*(rn-1),2*3*(rn-1))
    complex(kind=real64) :: Kl(2*3*(rn-1),2*3*(rn-1))
    ! Calculation of matrices
    L=fbem_fem_disrotra_L(rn,local_axis,nodal_axes)
    Kl=fbem_fem_disrotra_K_complex(rn,Kval)
    K=matmul(L,matmul(Kl,transpose(L)))
  end subroutine fbem_fem_disrotra_K_harmonic

  !! Calculation of all matrices needed for static analysis
  subroutine fbem_fem_disrotra_KLT_static(rn,local_axis,nodal_axes,Kval,KLT)
    implicit none
    ! I/O
    integer           :: rn                         !! Number of dimensions of the ambient space
    real(kind=real64) :: local_axis(rn,rn)          !! Local axes where springs and dashpots are defined
    real(kind=real64) :: nodal_axes(rn,rn,2)        !! Axes for each node DOFs nodal_axes(component,axis,node)
    real(kind=real64) :: Kval(8)                    !! Springs stiffnesses: 2D (k1,k2,kr3,k2r3), 3D (k1,k2,k3,kr1,kr2,kr3,k2r3,k3r2)
    real(kind=real64) :: KLT(2*3*(rn-1),2*3*(rn-1)) !! Stiffness matrix
    ! Local
    real(kind=real64) :: L(2*3*(rn-1),2*3*(rn-1))
    real(kind=real64) :: Kl(2*3*(rn-1),2*3*(rn-1))
    ! Calculation of matrices
    L=fbem_fem_disrotra_L(rn,local_axis,nodal_axes)
    Kl=fbem_fem_disrotra_K_real(rn,Kval)
    KLT=matmul(Kl,transpose(L))
  end subroutine fbem_fem_disrotra_KLT_static

  !! Calculation of all matrices needed for time harmonic analysis
  subroutine fbem_fem_disrotra_KLT_harmonic(rn,local_axis,nodal_axes,Kval,KLT)
    implicit none
    ! I/O
    integer              :: rn                         !! Number of dimensions of the ambient space
    real(kind=real64)    :: local_axis(rn,rn)          !! Local axes where springs and dashpots are defined
    real(kind=real64)    :: nodal_axes(rn,rn,2)        !! Axes for each node DOFs nodal_axes(component,axis,node)
    complex(kind=real64) :: Kval(8)                    !! Springs stiffnesses: 2D (k1,k2,kr3,k2r3), 3D (k1,k2,k3,kr1,kr2,kr3,k2r3,k3r2)
    complex(kind=real64) :: KLT(2*3*(rn-1),2*3*(rn-1)) !! Complex stiffness matrix
    ! Local
    real(kind=real64)    :: L(2*3*(rn-1),2*3*(rn-1))
    complex(kind=real64) :: Kl(2*3*(rn-1),2*3*(rn-1))
    ! Calculation of matrices
    L=fbem_fem_disrotra_L(rn,local_axis,nodal_axes)
    Kl=fbem_fem_disrotra_K_complex(rn,Kval)
    KLT=matmul(Kl,transpose(L))
  end subroutine fbem_fem_disrotra_KLT_harmonic

!~   ! Tests
!~   subroutine fbem_fem_disrotra_test()
!~     implicit none
!~     integer, parameter   :: rn=2
!~     real(kind=real64)    :: local_axes(rn,rn)   !! Local axes where springs and dashpots are defined
!~     real(kind=real64)    :: nodal_axes(rn,rn,2) !! Axes for each node DOFs nodal_axes(component,axis,node)
!~     real(kind=real64)    :: Kval(4*(rn-1))
!~     real(kind=real64)    :: K(2*3*(rn-1),2*3*(rn-1)) !! Complex stiffness matrix
!~     integer              :: kn
!~     real(kind=real64)    :: angle
!~     ! Static test
!~     angle=c_pi/2.d0
!~     local_axes(:,1)=[cos(angle),sin(angle)]!,0.d0]
!~     local_axes(:,2)=[-sin(angle),cos(angle)]!,0.d0]
!~     !local_axes(:,3)=[0.d0,0.d0,1.d0]
!~     nodal_axes(:,1,1)=local_axes(:,1)
!~     nodal_axes(:,2,1)=local_axes(:,2)
!~     !nodal_axes(:,3,1)=local_axes(:,3)
!~     nodal_axes(:,1,2)=local_axes(:,1)
!~     nodal_axes(:,2,2)=local_axes(:,2)
!~     !nodal_axes(:,3,2)=local_axes(:,3)
!~ !    nodal_axes(:,1,1)=[1,0,0]
!~ !    nodal_axes(:,2,1)=[0,1,0]
!~ !    nodal_axes(:,3,1)=[0,0,1]
!~ !    nodal_axes(:,1,2)=[1,0,0]
!~ !    nodal_axes(:,2,2)=[0,1,0]
!~ !    nodal_axes(:,3,2)=[0,0,1]
!~     nodal_axes(:,1,1)=[1,0]
!~     nodal_axes(:,2,1)=[0,1]
!~     nodal_axes(:,1,2)=[1,0]
!~     nodal_axes(:,2,2)=[0,1]
!~     Kval=[3,2,1,11]
!~     call fbem_fem_disrotra_K_static(rn,local_axes,nodal_axes,Kval,K)
!~     write(*,*) 'K'
!~     do kn=1,2*3*(rn-1)
!~       write(*,*) K(kn,:)
!~     end do
!~   end subroutine fbem_fem_disrotra_test

  ! ================================================================================================================================
  ! BAR ELEMENT
  ! ================================================================================================================================

  !! Coordinate transformation matrix L
  function fbem_fem_bar_L(rn,local_axis,nodal_axes)
    implicit none
    integer           :: rn                                    !! Number of dimensions of the ambient space
    real(kind=real64) :: local_axis(rn)                        !! Bar local axis
    real(kind=real64) :: nodal_axes(rn,rn,2)                   !! Axes for each node DOFs nodal_axes(component,axis,node)
    real(kind=real64) :: fbem_fem_bar_L(2*rn,2)
    integer           :: kn
    real(kind=real64) :: ep(rn,rn), e(rn,rn), Ln(rn,rn)
    fbem_fem_bar_L=0.d0
    ep=0.d0
    ep(:,1)=local_axis
    do kn=1,2
      e=nodal_axes(:,:,kn)
      call fbem_coordinate_transformation_L(rn,ep,e,Ln)
      fbem_fem_bar_L((rn*(kn-1)+1):(rn*kn),kn)=Ln(:,1)
    end do
  end function fbem_fem_bar_L

  !! Stiffness matrix for static analysis in local coordinates
  function fbem_fem_bar_K(L,A,E)
    implicit none
    real(kind=real64) :: L                   !! Length
    real(kind=real64) :: A                   !! Cross section area
    real(kind=real64) :: E                   !! Young' modulus
    real(kind=real64) :: fbem_fem_bar_K(2,2) !! Stiffness matrix in local coordinates
    fbem_fem_bar_K(1,1)= 1.d0
    fbem_fem_bar_K(1,2)=-1.d0
    fbem_fem_bar_K(2,1)=-1.d0
    fbem_fem_bar_K(2,2)= 1.d0
    fbem_fem_bar_K=fbem_fem_bar_K*E*A/L
  end function fbem_fem_bar_K

  !! Mass matrix in local coordinates
  function fbem_fem_bar_M(L,A,rho,Mtype)
    implicit none
    real(kind=real64) :: L                   !! Length
    real(kind=real64) :: A                   !! Cross section area
    real(kind=real64) :: rho                 !! Density
    integer           :: Mtype               !! Type of mass matrix: 0 (consistent), 1 (lumped mass)
    real(kind=real64) :: fbem_fem_bar_M(2,2) !! Mass matrix in local coordinates
    select case (Mtype)
      case (0)
        fbem_fem_bar_M(1,1)=1.0d0
        fbem_fem_bar_M(1,2)=0.5d0
        fbem_fem_bar_M(2,1)=0.5d0
        fbem_fem_bar_M(2,2)=1.0d0
        fbem_fem_bar_M=fbem_fem_bar_M*rho*A*L/3.d0
      case (1)
        fbem_fem_bar_M=0.d0
        fbem_fem_bar_M(1,1)=1.0d0
        fbem_fem_bar_M(2,2)=1.0d0
        fbem_fem_bar_M=fbem_fem_bar_M*rho*A*L/2.d0
      case default
        call fbem_error_message(error_unit,0,__FILE__,__LINE__,'Mtype should be 0 (consistent) or 1 (lumped)')
    end select
  end function fbem_fem_bar_M

  !! Calculation of all matrices needed for static analysis
  subroutine fbem_fem_bar_K_static(rn,x,A,nodal_axes,E,K)
    implicit none
    ! I/O
    integer           :: rn                  !! Number of dimensions of the ambient space
    real(kind=real64) :: x(rn,2)             !! Coordinates of the nodes
    real(kind=real64) :: A                   !! Cross section area
    real(kind=real64) :: nodal_axes(rn,rn,2) !! Axes for each node DOFs nodal_axes(component,axis,node)
    real(kind=real64) :: E                   !! Young' modulus
    real(kind=real64) :: K(2*rn,2*rn)        !! Stiffness matrix
    ! Local
    real(kind=real64) :: L(2*rn,2)
    real(kind=real64) :: Le, local_axis(rn)
    real(kind=real64) :: Kl(2,2)
    ! Auxiliary calculations
    local_axis=x(:,2)-x(:,1)
    Le=sqrt(dot_product(local_axis,local_axis))
    local_axis=local_axis/Le
    ! Calculation of matrices
    L=fbem_fem_bar_L(rn,local_axis,nodal_axes)
    Kl=fbem_fem_bar_K(Le,A,E)
    K=matmul(L,matmul(Kl,transpose(L)))
  end subroutine fbem_fem_bar_K_static

  !! Calculation of all matrices needed for time harmonic analysis
  subroutine fbem_fem_bar_K_harmonic(rn,omega,x,A,nodal_axes,E,rho,Mtype,K)
    implicit none
    ! I/O
    integer              :: rn                  !! Number of dimensions of the ambient space
    real(kind=real64)    :: omega               !! Circular frequency
    real(kind=real64)    :: x(rn,2)             !! Coordinates of the nodes
    real(kind=real64)    :: A                   !! Cross section area
    real(kind=real64)    :: nodal_axes(rn,rn,2) !! Axes for each node DOFs nodal_axes(component,axis,node)
    complex(kind=real64) :: E                   !! Young' modulus
    real(kind=real64)    :: rho                 !! Density
    integer              :: Mtype               !! Type of mass matrix: 0 (consistent), 1 (lumped)
    complex(kind=real64) :: K(2*rn,2*rn)        !! Stiffness matrix
    ! Local
    real(kind=real64)    :: L(2*rn,2)
    real(kind=real64)    :: Le, local_axis(rn)
    complex(kind=real64) :: Kl(2,2)
    real(kind=real64)    :: Ml(2,2)
    ! Auxiliary calculations
    local_axis=x(:,2)-x(:,1)
    Le=sqrt(dot_product(local_axis,local_axis))
    local_axis=local_axis/Le
    ! Calculation of matrices
    L=fbem_fem_bar_L(rn,local_axis,nodal_axes)
    Kl=fbem_fem_bar_K(Le,A,1.d0)
    Kl=E*Kl
    Ml=fbem_fem_bar_M(Le,A,rho,Mtype)
    Kl=Kl-omega**2*Ml
    K=matmul(L,matmul(Kl,transpose(L)))
  end subroutine fbem_fem_bar_K_harmonic

 !! D*B product for stress resultant calculation
  function fbem_fem_bar_DB(etype,L,A,E,xi)
    implicit none
    integer           :: etype               !! Type of element
    real(kind=real64) :: L                   !! Length
    real(kind=real64) :: A                   !! Cross section area
    real(kind=real64) :: E                   !! oung' modulus
    real(kind=real64) :: xi                  !! Reference coordinate (0<=xi<=1)
    real(kind=real64) :: fbem_fem_bar_DB(fbem_n_nodes(etype))
    fbem_fem_bar_DB=0
    select case (etype)
      case (fbem_line2)
        fbem_fem_bar_DB(1)=-1
        fbem_fem_bar_DB(2)= 1
      case (fbem_line3)
        fbem_fem_bar_DB(1)=-3+4*xi
        fbem_fem_bar_DB(2)=-1+4*xi
        fbem_fem_bar_DB(3)=4-8*xi
      case default
        call fbem_error_message(error_unit,0,__FILE__,__LINE__,'invalid type of bar')
    end select
    fbem_fem_bar_DB=E*A/L*fbem_fem_bar_DB
  end function fbem_fem_bar_DB

  ! Stress resultants extrapolated from Gauss points. Remember Fsigma must be later multiplied by E (Young modulus)
  subroutine fbem_fem_bar_stress_resultants(rn,x,A,nodal_axes,setype,sedelta,Fsigma)
    implicit none
    ! I/O
    integer           :: rn                          !! Number of dimensions of the ambient space
    real(kind=real64) :: x(rn,2)                     !! Coordinates of the nodes
    real(kind=real64) :: A                           !! Length
    real(kind=real64) :: nodal_axes(rn,rn,2)         !! Axes for each node DOFs nodal_axes(component,axis,node)
    integer           :: setype                      !! Stress interpolation: type of interpolation
    real(kind=real64) :: sedelta                     !! Stress interpolation: delta of type of interpolation
    real(kind=real64), allocatable :: Fsigma(:,:,:)  !! Stress resultants matrix at interpolation points (without multiplying by E)
    ! Local
    integer           :: ksp
    real(kind=real64) :: Le, local_axis(rn)
    real(kind=real64) :: L(2*rn,2), DB(1,2)
    real(kind=real64) :: xi
    !
    ! Initialization
    !
    local_axis=x(:,2)-x(:,1)
    Le=sqrt(dot_product(local_axis,local_axis))
    local_axis=local_axis/Le
    ! Coordinate transformation matrix
    L=fbem_fem_bar_L(rn,local_axis,nodal_axes)
    ! Select the stress interpolation scheme
    setype=fbem_line1
    sedelta=0
    ! Local stress resultants matrix at interpolation points
    allocate(Fsigma(1,2*rn,fbem_n_nodes(setype)))
    Fsigma=0
    !
    ! Loops through sampling points
    !
    do ksp=1,fbem_n_nodes(setype)
      !
      ! Sampling point coordinate
      !
#     define node ksp
#     define etype setype
#     define delta sedelta
#     include <xi_1d_at_node.rc>
#     undef node
#     undef etype
#     undef delta
      !
      ! Convert -1<=xi<=1 to 0<=xi<=1
      !
      DB(1,:)=fbem_fem_bar_DB(fbem_line2,Le,A,1.d0,xi)
      Fsigma(:,:,ksp)=matmul(DB,transpose(L))
    end do
  end subroutine fbem_fem_bar_stress_resultants

!~   ! Tests
!~   subroutine fbem_fem_bar_test()
!~     implicit none
!~     integer, parameter   :: rn=2
!~     real(kind=real64)    :: x(rn,2)             !! Coordinates of the nodes
!~     real(kind=real64)    :: nodal_axes(rn,rn,2) !! Axes for each node DOFs nodal_axes(component,axis,node)
!~     real(kind=real64)    :: K(2*rn,2*rn)        !! Complex stiffness matrix
!~     integer              :: kn
!~     complex(kind=real64) :: Kc(2*rn,2*rn)       !! Stiffness matrix
!~     ! Static test
!~     x(:,1)=[0,0]
!~     x(:,2)=[1,0]
!~     nodal_axes(:,1,1)=[1,0]
!~     nodal_axes(:,2,1)=[0,1]
!~     nodal_axes(:,1,2)=[1,0]
!~     nodal_axes(:,2,2)=[0,1]
!~ !    nodal_axes(:,1,1)=[1,0,0]
!~ !    nodal_axes(:,2,1)=[0,1,0]
!~ !    nodal_axes(:,3,1)=[0,0,1]
!~ !    nodal_axes(:,1,2)=[1,0,0]
!~ !    nodal_axes(:,2,2)=[0,1,0]
!~ !    nodal_axes(:,3,2)=[0,0,1]
!~     call fbem_fem_bar_K_static(rn,x,1.d0,nodal_axes,1.d0,K)
!~     write(*,*) 'K'
!~     do kn=1,2*rn
!~       write(*,*) K(kn,:)
!~     end do
!~     call fbem_fem_bar_K_harmonic(rn,2.d0*sqrt(3.d0),x,1.d0,nodal_axes,(1.d0,0.d0),1.d0,1,Kc)
!~     write(*,*) 'Kc'
!~     do kn=1,2*rn
!~       write(*,*) dreal(Kc(kn,:))
!~     end do
!~     call fbem_fem_bar_K_harmonic(rn,2.d0,x,1.d0,nodal_axes,(1.d0,0.d0),1.d0,2,Kc)
!~     write(*,*) 'Kc'
!~     do kn=1,2*rn
!~       write(*,*) dreal(Kc(kn,:))
!~     end do
!~   end subroutine fbem_fem_bar_test

  ! ================================================================================================================================
  ! STRAIGHT EULER-BERNOULLI AND TIMOSHENKO BEAM ELEMENTS
  ! ================================================================================================================================

  !! Coordinate transformation matrix L for all straight beam elements, from local to an specified coordinate system for each node.
  !! Calculation of transformation matrix of a=L*a'
  function fbem_fem_strbeam_L_element(rn,etype,local_axis,nodal_axes)
    implicit none
    integer           :: rn                                    !! Number of dimensions of the ambient space
    integer           :: etype                                 !! Type of element
    real(kind=real64) :: local_axis(rn,rn)                     !! Beam local axis
    real(kind=real64) :: nodal_axes(rn,rn,fbem_n_nodes(etype)) !! Axes for each node DOFs nodal_axes(component,axis,node)
    real(kind=real64) :: fbem_fem_strbeam_L_element(3*(rn-1)*fbem_n_nodes(etype),3*(rn-1)*fbem_n_nodes(etype))
    integer           :: kn, ki, kj
    real(kind=real64) :: e(rn,rn), Ln(rn,rn)
    fbem_fem_strbeam_L_element=0
    select case (etype)
      case (fbem_line2)
      case (fbem_line3)
      case default
        call fbem_error_message(error_unit,0,__FILE__,__LINE__,'invalid type of strbeam')
    end select
    ki=1
    do kn=1,fbem_n_nodes(etype)
      e=nodal_axes(:,:,kn)
      call fbem_coordinate_transformation_L(rn,local_axis,e,Ln)
      kj=ki+rn-1
      fbem_fem_strbeam_L_element(ki:kj,ki:kj)=Ln
      if (rn.eq.2) then
        fbem_fem_strbeam_L_element(kj+1,kj+1)=1
      else
        fbem_fem_strbeam_L_element((kj+1):(kj+3),(kj+1):(kj+3))=Ln
      end if
      ki=ki+3*(rn-1)
    end do
  end function fbem_fem_strbeam_L_element

  !! Calculate the shape functions matrix N at a given local coordinate
  function fbem_fem_strbeam_Na(etype,xi) result (N)
    implicit none
    ! I/O
    integer           :: etype                     !! Type of element (displacements interpolation)
    real(kind=real64) :: xi                        !! Local coordinate
    real(kind=real64) :: N(fbem_n_nodes(etype))    !! Shape functions matrix (DOF in local coordinates)
    ! Local
    real(kind=real64) :: aux(10)                   ! Auxiliary variable needed for shape_functions module resources
    N=0
#   define phi N
#   define delta 0.d0
#   include <phi_1d.rc>
#   undef delta
#   undef phi
  end function fbem_fem_strbeam_Na

  !! Axial stiffness matrix (displacement or torsion)
  function fbem_fem_strbeam_Ka(etype,L,A,E)
    implicit none
    integer           :: etype               !! Type of element
    real(kind=real64) :: L                   !! Length
    real(kind=real64) :: A                   !! For axial displacement: Cross section area; For torsion: torsional constant.
    real(kind=real64) :: E                   !! For axial displacement: Young' modulus; For torsion: shear modulus
    real(kind=real64) :: fbem_fem_strbeam_Ka(fbem_n_nodes(etype),fbem_n_nodes(etype))
    fbem_fem_strbeam_Ka=0
    select case (etype)
      case (fbem_line2)
        fbem_fem_strbeam_Ka(1,1)= 1
        fbem_fem_strbeam_Ka(2,1)=-1
        fbem_fem_strbeam_Ka(1,2)=-1
        fbem_fem_strbeam_Ka(2,2)= 1
        fbem_fem_strbeam_Ka=fbem_fem_strbeam_Ka*E*A/L
      case (fbem_line3)
        fbem_fem_strbeam_Ka(1,1)= 7
        fbem_fem_strbeam_Ka(2,1)= 1
        fbem_fem_strbeam_Ka(3,1)=-8
        fbem_fem_strbeam_Ka(1,2)= 1
        fbem_fem_strbeam_Ka(2,2)= 7
        fbem_fem_strbeam_Ka(3,2)=-8
        fbem_fem_strbeam_Ka(1,3)=-8
        fbem_fem_strbeam_Ka(2,3)=-8
        fbem_fem_strbeam_Ka(3,3)=16
        fbem_fem_strbeam_Ka=fbem_fem_strbeam_Ka*E*A/(3*L)
      case default
        call fbem_error_message(error_unit,0,__FILE__,__LINE__,'invalid type of strbeam')
    end select
  end function fbem_fem_strbeam_Ka

  !! Axial mass matrix (displacement or torsion)
  function fbem_fem_strbeam_Ma(etype,L,A,rho)
    implicit none
    integer           :: etype               !! Type of element
    real(kind=real64) :: L                   !! Length
    real(kind=real64) :: A                   !! For axial displacement: Cross section area; For torsion: torsion moment of inertia
    real(kind=real64) :: rho                 !! Material density
    real(kind=real64) :: fbem_fem_strbeam_Ma(fbem_n_nodes(etype),fbem_n_nodes(etype))
    fbem_fem_strbeam_Ma=0
    select case (etype)
      case (fbem_line2)
        fbem_fem_strbeam_Ma(1,1)=2
        fbem_fem_strbeam_Ma(2,1)=1
        fbem_fem_strbeam_Ma(1,2)=1
        fbem_fem_strbeam_Ma(2,2)=2
        fbem_fem_strbeam_Ma=fbem_fem_strbeam_Ma*rho*A*L/6
      case (fbem_line3)
        fbem_fem_strbeam_Ma(1,1)= 4
        fbem_fem_strbeam_Ma(2,1)=-1
        fbem_fem_strbeam_Ma(3,1)= 2
        fbem_fem_strbeam_Ma(1,2)=-1
        fbem_fem_strbeam_Ma(2,2)= 4
        fbem_fem_strbeam_Ma(3,2)= 2
        fbem_fem_strbeam_Ma(1,3)= 2
        fbem_fem_strbeam_Ma(2,3)= 2
        fbem_fem_strbeam_Ma(3,3)=16
        fbem_fem_strbeam_Ma=fbem_fem_strbeam_Ma*rho*A*L/30
      case default
        call fbem_error_message(error_unit,0,__FILE__,__LINE__,'invalid type of strbeam')
    end select
  end function fbem_fem_strbeam_Ma

  !! Axial equivalent forces matrix (displacement or torsion)
  function fbem_fem_strbeam_Qa(etype,L)
    implicit none
    integer           :: etype               !! Type of element
    real(kind=real64) :: L                   !! Length
    real(kind=real64) :: fbem_fem_strbeam_Qa(fbem_n_nodes(etype),fbem_n_nodes(etype))
    fbem_fem_strbeam_Qa=0
    select case (etype)
      case (fbem_line2)
        fbem_fem_strbeam_Qa(1,1)=2
        fbem_fem_strbeam_Qa(2,1)=1
        fbem_fem_strbeam_Qa(1,2)=1
        fbem_fem_strbeam_Qa(2,2)=2
        fbem_fem_strbeam_Qa=fbem_fem_strbeam_Qa*L/6
      case (fbem_line3)
        fbem_fem_strbeam_Qa(1,1)= 4
        fbem_fem_strbeam_Qa(2,1)=-1
        fbem_fem_strbeam_Qa(3,1)= 2
        fbem_fem_strbeam_Qa(1,2)=-1
        fbem_fem_strbeam_Qa(2,2)= 4
        fbem_fem_strbeam_Qa(3,2)= 2
        fbem_fem_strbeam_Qa(1,3)= 2
        fbem_fem_strbeam_Qa(2,3)= 2
        fbem_fem_strbeam_Qa(3,3)=16
        fbem_fem_strbeam_Qa=fbem_fem_strbeam_Qa*L/30
      case default
        call fbem_error_message(error_unit,0,__FILE__,__LINE__,'invalid type of strbeam')
    end select
  end function fbem_fem_strbeam_Qa

  !! Calculate the shape functions matrix N at a given local coordinate
  function fbem_fem_strbeam_Nl(etype,esubtype,L,phi,xi) result (N)
    implicit none
    ! I/O
    integer           :: etype                     !! Type of element (displacements interpolation)
    integer           :: esubtype                  !! Subtype (for line3)
    real(kind=real64) :: L                         !! Length
    real(kind=real64) :: phi                       !! Ratio of the beam bending stiffness to the shear stiffness
    real(kind=real64) :: xi                        !! Local coordinate (between -1 and 1)
    real(kind=real64) :: N(2*fbem_n_nodes(etype))  !! Shape functions matrix (DOF in local coordinates)
    ! Local
    real(kind=real64) :: xip                       !! Local coordinate (between 0 and 1)
    N=0
    xip=0.5d0*(xi+1.d0)
    select case (etype)
      case (fbem_line2)
        N(1)=-((-1.d0+xip)*(1.d0+xip-2.d0*xip**2+phi))/(1.d0+phi)
        N(2)=-(L*(-1.d0+xip)*xip*(2.d0-2.d0*xip+phi))/(2.d0*(1.d0+phi))
        N(3)=(xip*((3.d0-2.d0*xip)*xip+phi))/(1.d0+phi)
        N(4)=(L*(-1.d0+xip)*xip*(2.d0*xip+phi))/(2.d0*(1.d0+phi))
      case (fbem_line3)
        if (esubtype.eq.0) then
          N(1)=-(((-1.d0+xip)*(-1.d0+2.d0*xip)*(4.d0*xip**2*(1.d0+phi)-3.d0*xip*(1.d0+2.d0*phi)-(1.d0+phi)*(1.d0+4.d0*phi)))/((1.d0+phi)*(1.d0+4.d0*phi)))
          N(2)=-(L*xip*(1.d0-3.d0*xip+2.d0*xip**2)*(-2.d0-5.d0*phi+2.d0*xip*(1.d0+phi)))/(2.d0+10.d0*phi+8.d0*phi**2)
          N(3)=-((xip*(-1.d0+2.d0*xip)*(4.d0*xip**2*(1.d0+phi)-xip*(5.d0+2.d0*phi)-phi*(7.d0+4.d0*phi)))/(1.d0+5.d0*phi+4.d0*phi**2))
          N(4)=((L*xip*(1.d0-3.d0*xip+2.d0*xip**2)*(3.d0*phi+2.d0*xip*(1.d0+phi)))/(2.d0+10.d0*phi+8.d0*phi**2))
          N(5)=(16.d0*(-1.d0+xip)*xip*((-1.d0+xip)*xip-phi))/(1.d0+4.d0*phi)
        else
          N(1)=((-1.d0+xip)*(-1.d0+2.d0*xip)*(1.d0+xip*(3.d0+4.d0*xip*(-4.d0+3.d0*xip))+9.d0*phi+2.d0*xip*(-3.d0+2.d0*xip)*(1.d0+12.d0*xip)*phi+20.d0*(1.d0-4.d0*xip)*phi**2))/((1.d0+4.d0*phi)*(1.d0+5.d0*phi))
          N(2)=-(L*(-1.d0+xip)*xip*(-1.d0+2.d0*xip)*(-6.d0+12.d0*xip**2*(-1.d0+2.d0*phi)*(1.d0+4.d0*phi)+phi*(-15.d0+8.d0*(1.d0-20.d0*phi)*phi)+6.d0*xip*(3.d0+(9.d0-16.d0*phi)*phi)))/(6.d0*(1.d0+4.d0*phi)*(1.d0+5.d0*phi))
          N(3)=-((xip*(-1.d0+2.d0*xip)*(12.d0*xip**3*(1.d0+4.d0*phi)-4.d0*xip**2*(5.d0+19.d0*phi)+phi*(17.d0+60.d0*phi)+xip*(7.d0+2.d0*phi-80.d0*phi**2)))/(1.d0+9.d0*phi+20.d0*phi**2))
          N(4)=-(L*(-1.d0+xip)*xip*(-1.d0+2.d0*xip)*(12.d0*xip**2*(-1.d0+2.d0*phi)*(1.d0+4.d0*phi)+phi*(15.d0+8.d0*(1.d0-20.d0*phi)*phi)-6.d0*xip*(-1.d0+phi+16.d0*phi**2)))/(6.d0*(1.d0+4.d0*phi)*(1.d0+5.d0*phi))
          N(5)=(16.d0*(-1.d0+xip)*xip*((-1.d0+xip)*xip-phi))/(1.d0+4.d0*phi)
          N(6)=-(-8.d0*L*(-1.d0+xip)*xip*(-1.d0+2.d0*xip)*(3.d0*(-1.d0+xip)*xip+3.d0*(-2.d0+xip)*(1.d0+xip)*phi-5.d0*phi**2))/(3.d0+15.d0*phi)
        end if
      case default
        call fbem_error_message(error_unit,0,__FILE__,__LINE__,'invalid type of strbeam')
    end select
  end function fbem_fem_strbeam_Nl

  !! Lateral stiffness matrix
  function fbem_fem_strbeam_Kl(etype,esubtype,L,I,E,phi)
    implicit none
    integer           :: etype               !! Type of element
    integer           :: esubtype            !! Subtype (for line3)
    real(kind=real64) :: L                   !! Length
    real(kind=real64) :: I                   !! Moment of inertia
    real(kind=real64) :: E                   !! Young' modulus
    real(kind=real64) :: phi                 !! Ratio of the beam bending stiffness to the shear stiffness
    real(kind=real64) :: fbem_fem_strbeam_Kl(2*fbem_n_nodes(etype),2*fbem_n_nodes(etype))
    fbem_fem_strbeam_Kl=0
    select case (etype)
      case (fbem_line2)
        fbem_fem_strbeam_Kl(1,1)=12
        fbem_fem_strbeam_Kl(2,1)=6*L
        fbem_fem_strbeam_Kl(3,1)=-12
        fbem_fem_strbeam_Kl(4,1)=6*L
        fbem_fem_strbeam_Kl(1,2)=6*L
        fbem_fem_strbeam_Kl(2,2)=L**2*(4+phi)
        fbem_fem_strbeam_Kl(3,2)=-6*L
        fbem_fem_strbeam_Kl(4,2)=-(L**2*(-2+phi))
        fbem_fem_strbeam_Kl(1,3)=-12
        fbem_fem_strbeam_Kl(2,3)=-6*L
        fbem_fem_strbeam_Kl(3,3)=12
        fbem_fem_strbeam_Kl(4,3)=-6*L
        fbem_fem_strbeam_Kl(1,4)=6*L
        fbem_fem_strbeam_Kl(2,4)=-(L**2*(-2+phi))
        fbem_fem_strbeam_Kl(3,4)=-6*L
        fbem_fem_strbeam_Kl(4,4)=L**2*(4+phi)
        fbem_fem_strbeam_Kl=fbem_fem_strbeam_Kl*E*I/(1+phi)/L**3
      case (fbem_line3)
        if (esubtype.eq.0) then
          fbem_fem_strbeam_Kl(1,1)=4*(79+56*phi*(9+10*phi))
          fbem_fem_strbeam_Kl(2,1)=2*L*(47+8*phi*(39+50*phi))
          fbem_fem_strbeam_Kl(3,1)=4*(49+8*phi*(33+10*phi))
          fbem_fem_strbeam_Kl(4,1)=2*L*(-17+8*phi*(-9+10*phi))
          fbem_fem_strbeam_Kl(5,1)=-512*(1+phi)*(1+5*phi)
          fbem_fem_strbeam_Kl(1,2)=2*L*(47+8*phi*(39+50*phi))
          fbem_fem_strbeam_Kl(2,2)=L**2*(36+phi*(261+40*phi*(11+2*phi)))
          fbem_fem_strbeam_Kl(3,2)=2*L*(17+8*(9-10*phi)*phi)
          fbem_fem_strbeam_Kl(4,2)=-(L**2*(6+phi*(21+40*phi*(-1+2*phi))))
          fbem_fem_strbeam_Kl(5,2)=-128*L*(1+phi)*(1+5*phi)
          fbem_fem_strbeam_Kl(1,3)=4*(49+8*phi*(33+10*phi))
          fbem_fem_strbeam_Kl(2,3)=2*L*(17+8*(9-10*phi)*phi)
          fbem_fem_strbeam_Kl(3,3)=4*(79+56*phi*(9+10*phi))
          fbem_fem_strbeam_Kl(4,3)=-2*L*(47+8*phi*(39+50*phi))
          fbem_fem_strbeam_Kl(5,3)=-512*(1+phi)*(1+5*phi)
          fbem_fem_strbeam_Kl(1,4)=2*L*(-17+8*phi*(-9+10*phi))
          fbem_fem_strbeam_Kl(2,4)=-(L**2*(6+phi*(21+40*phi*(-1+2*phi))))
          fbem_fem_strbeam_Kl(3,4)=-2*L*(47+8*phi*(39+50*phi))
          fbem_fem_strbeam_Kl(4,4)=L**2*(36+phi*(261+40*phi*(11+2*phi)))
          fbem_fem_strbeam_Kl(5,4)=128*L*(1+phi)*(1+5*phi)
          fbem_fem_strbeam_Kl(1,5)=-512*(1+phi)*(1+5*phi)
          fbem_fem_strbeam_Kl(2,5)=-128*L*(1+phi)*(1+5*phi)
          fbem_fem_strbeam_Kl(3,5)=-512*(1+phi)*(1+5*phi)
          fbem_fem_strbeam_Kl(4,5)=128*L*(1+phi)*(1+5*phi)
          fbem_fem_strbeam_Kl(5,5)=1024*(1+phi)*(1+5*phi)
          fbem_fem_strbeam_Kl=fbem_fem_strbeam_Kl*E*I/5/(1+phi)/(1+4*phi)**2/L**3
        else
          fbem_fem_strbeam_Kl(1,1)=12*(1273+5*phi*(4869+40*phi*(675+1162*phi)))
          fbem_fem_strbeam_Kl(2,1)=6*L*(569+5*phi*(1989+8*phi*(1183+10*(111-224*phi)*phi)))
          fbem_fem_strbeam_Kl(3,1)=-12*(377+5*phi*(2181+40*phi*(339+602*phi)))
          fbem_fem_strbeam_Kl(4,1)=6*L*(121+5*phi*(645-8*phi*(-343+10*phi*(29+224*phi))))
          fbem_fem_strbeam_Kl(5,1)=-10752*(1+5*phi)**3
          fbem_fem_strbeam_Kl(6,1)=1920*L*(1+4*phi)**2*(3+5*phi*(9+7*phi))
          fbem_fem_strbeam_Kl(1,2)=6*L*(569+5*phi*(1989+8*phi*(1183+10*(111-224*phi)*phi)))
          fbem_fem_strbeam_Kl(2,2)=L**2*(996+5*phi*(3321+phi*(16133+8*phi*(2769+10*phi*(9+448*phi)))))
          fbem_fem_strbeam_Kl(3,2)=6*L*(-121+5*phi*(-645+8*phi*(-343+10*phi*(29+224*phi))))
          fbem_fem_strbeam_Kl(4,2)=L**2*(114+5*phi*(549+phi*(971+8*phi*(-1221+10*phi*(-201+448*phi)))))
          fbem_fem_strbeam_Kl(5,2)=-2688*L*(1+5*phi)**3
          fbem_fem_strbeam_Kl(6,2)=-320*(-1+2*phi)*(L+4*L*phi)**2*(3+5*phi*(9+7*phi))
          fbem_fem_strbeam_Kl(1,3)=-12*(377+5*phi*(2181+40*phi*(339+602*phi)))
          fbem_fem_strbeam_Kl(2,3)=6*L*(-121+5*phi*(-645+8*phi*(-343+10*phi*(29+224*phi))))
          fbem_fem_strbeam_Kl(3,3)=12*(1273+5*phi*(4869+40*phi*(675+1162*phi)))
          fbem_fem_strbeam_Kl(4,3)=6*L*(-569+5*phi*(-1989+8*phi*(-1183+10*phi*(-111+224*phi))))
          fbem_fem_strbeam_Kl(5,3)=-10752*(1+5*phi)**3
          fbem_fem_strbeam_Kl(6,3)=-1920*L*(1+4*phi)**2*(3+5*phi*(9+7*phi))
          fbem_fem_strbeam_Kl(1,4)=6*L*(121+5*phi*(645-8*phi*(-343+10*phi*(29+224*phi))))
          fbem_fem_strbeam_Kl(2,4)=L**2*(114+5*phi*(549+phi*(971+8*phi*(-1221+10*phi*(-201+448*phi)))))
          fbem_fem_strbeam_Kl(3,4)=6*L*(-569+5*phi*(-1989+8*phi*(-1183+10*phi*(-111+224*phi))))
          fbem_fem_strbeam_Kl(4,4)=L**2*(996+5*phi*(3321+phi*(16133+8*phi*(2769+10*phi*(9+448*phi)))))
          fbem_fem_strbeam_Kl(5,4)=2688*L*(1+5*phi)**3
          fbem_fem_strbeam_Kl(6,4)=-320*(-1+2*phi)*(L+4*L*phi)**2*(3+5*phi*(9+7*phi))
          fbem_fem_strbeam_Kl(1,5)=-10752*(1+5*phi)**3
          fbem_fem_strbeam_Kl(2,5)=-2688*L*(1+5*phi)**3
          fbem_fem_strbeam_Kl(3,5)=-10752*(1+5*phi)**3
          fbem_fem_strbeam_Kl(4,5)=2688*L*(1+5*phi)**3
          fbem_fem_strbeam_Kl(5,5)=21504*(1+5*phi)**3
          fbem_fem_strbeam_Kl(6,5)=0
          fbem_fem_strbeam_Kl(1,6)=1920*L*(1+4*phi)**2*(3+5*phi*(9+7*phi))
          fbem_fem_strbeam_Kl(2,6)=-320*(-1+2*phi)*(L+4*L*phi)**2*(3+5*phi*(9+7*phi))
          fbem_fem_strbeam_Kl(3,6)=-1920*L*(1+4*phi)**2*(3+5*phi*(9+7*phi))
          fbem_fem_strbeam_Kl(4,6)=-320*(-1+2*phi)*(L+4*L*phi)**2*(3+5*phi*(9+7*phi))
          fbem_fem_strbeam_Kl(5,6)=0
          fbem_fem_strbeam_Kl(6,6)=1280*(1+phi)*(L+4*L*phi)**2*(3+5*phi*(9+7*phi))
          fbem_fem_strbeam_Kl=fbem_fem_strbeam_Kl*E*I/105/(1+5*phi)**2/(1+4*phi)**2/L**3
        end if
      case default
        call fbem_error_message(error_unit,0,__FILE__,__LINE__,'invalid type of strbeam')
    end select
  end function fbem_fem_strbeam_Kl

  !! Lateral mass matrix
  function fbem_fem_strbeam_Ml(etype,esubtype,L,A,I,rho,phi)
    implicit none
    integer           :: etype               !! Type of element
    integer           :: esubtype            !! Subtype (for line3)
    real(kind=real64) :: L                   !! Length
    real(kind=real64) :: A                   !! Cross-section area
    real(kind=real64) :: I                   !! Moment of inertia
    real(kind=real64) :: rho                 !! Density
    real(kind=real64) :: phi                 !! Ratio of the beam bending stiffness to the shear stiffness
    real(kind=real64) :: fbem_fem_strbeam_Ml(2*fbem_n_nodes(etype),2*fbem_n_nodes(etype))
    real(kind=real64) :: Mt(2*fbem_n_nodes(etype),2*fbem_n_nodes(etype))
    real(kind=real64) :: Mr(2*fbem_n_nodes(etype),2*fbem_n_nodes(etype))
    Mr=0
    Mt=0
    select case (etype)
      case (fbem_line2)
        Mt(1,1)=4*(78+7*phi*(21+10*phi))
        Mt(2,1)=L*(44+7*phi*(11+5*phi))
        Mt(3,1)=4*(27+7*phi*(9+5*phi))
        Mt(4,1)=-(L*(26+7*phi*(9+5*phi)))
        Mt(1,2)=L*(44+7*phi*(11+5*phi))
        Mt(2,2)=L**2*(8+7*phi*(2+phi))
        Mt(3,2)=L*(26+7*phi*(9+5*phi))
        Mt(4,2)=-(L**2*(6+7*phi*(2+phi)))
        Mt(1,3)=4*(27+7*phi*(9+5*phi))
        Mt(2,3)=L*(26+7*phi*(9+5*phi))
        Mt(3,3)=4*(78+7*phi*(21+10*phi))
        Mt(4,3)=-(L*(44+7*phi*(11+5*phi)))
        Mt(1,4)=-(L*(26+7*phi*(9+5*phi)))
        Mt(2,4)=-(L**2*(6+7*phi*(2+phi)))
        Mt(3,4)=-(L*(44+7*phi*(11+5*phi)))
        Mt(4,4)=L**2*(8+7*phi*(2+phi))
        Mt=Mt*rho*A*L/(1+phi)**2/840
        Mr(1,1)=36
        Mr(2,1)=-3*L*(-1+5*phi)
        Mr(3,1)=-36
        Mr(4,1)=-3*L*(-1+5*phi)
        Mr(1,2)=-3*L*(-1+5*phi)
        Mr(2,2)=L**2*(4+5*phi*(1+2*phi))
        Mr(3,2)=-3*L*(1-5*phi)
        Mr(4,2)=L**2*(-1+5*(-1+phi)*phi)
        Mr(1,3)=-36
        Mr(2,3)=-3*L*(1-5*phi)
        Mr(3,3)=36
        Mr(4,3)=-3*L*(1-5*phi)
        Mr(1,4)=-3*L*(-1+5*phi)
        Mr(2,4)=L**2*(-1+5*(-1+phi)*phi)
        Mr(3,4)=-3*L*(1-5*phi)
        Mr(4,4)=L**2*(4+5*phi*(1+2*phi))
        Mr=Mr*rho*I/(1+phi)**2/30/L
      case (fbem_line3)
        if (esubtype.eq.0) then
          Mt(1,1)=4*(130+phi*(1169+2*phi*(1745+12*phi*(157+56*phi))))
          Mt(2,1)=L*(40+phi*(299+phi*(685+372*phi)))
          Mt(3,1)=-4*(23+phi*(307+phi*(1079+48*phi*(23+7*phi))))
          Mt(4,1)=L*(14+5*phi*(5+4*phi)*(7+15*phi))
          Mt(5,1)=32*(1+phi)**2*(5+45*phi+84*phi**2)
          Mt(1,2)=L*(40+phi*(299+phi*(685+372*phi)))
          Mt(2,2)=L**2*(4+phi*(26+49*phi))
          Mt(3,2)=-(L*(14+5*phi*(5+4*phi)*(7+15*phi)))
          Mt(4,2)=L**2*(2+phi*(22+47*phi))
          Mt(5,2)=16*L*(1+phi)**2*(1+6*phi)
          Mt(1,3)=-4*(23+phi*(307+phi*(1079+48*phi*(23+7*phi))))
          Mt(2,3)=-(L*(14+5*phi*(5+4*phi)*(7+15*phi)))
          Mt(3,3)=4*(130+phi*(1169+2*phi*(1745+12*phi*(157+56*phi))))
          Mt(4,3)=-(L*(40+phi*(299+phi*(685+372*phi))))
          Mt(5,3)=32*(1+phi)**2*(5+45*phi+84*phi**2)
          Mt(1,4)=L*(14+5*phi*(5+4*phi)*(7+15*phi))
          Mt(2,4)=L**2*(2+phi*(22+47*phi))
          Mt(3,4)=-(L*(40+phi*(299+phi*(685+372*phi))))
          Mt(4,4)=L**2*(4+phi*(26+49*phi))
          Mt(5,4)=-16*L*(1+phi)**2*(1+6*phi)
          Mt(1,5)=32*(1+phi)**2*(5+45*phi+84*phi**2)
          Mt(2,5)=16*L*(1+phi)**2*(1+6*phi)
          Mt(3,5)=32*(1+phi)**2*(5+45*phi+84*phi**2)
          Mt(4,5)=-16*L*(1+phi)**2*(1+6*phi)
          Mt(5,5)=1024*(1+phi)**2*(1+3*phi*(3+7*phi))
          Mt=Mt*rho*A*L/(1+phi)**2/(1+4*phi)**2/2520
          Mr(1,1)=508+32*phi*(79+134*phi)
          Mr(2,1)=-(L*(-29+phi*(145+16*phi*(59+119*phi))))
          Mr(3,1)=-4*(-1+8*phi*(47+118*phi))
          Mr(4,1)=-(L*(-13+phi*(-271+16*phi*(4+91*phi))))
          Mr(5,1)=-512*(1+phi)**2
          Mr(1,2)=-(L*(-29+phi*(145+16*phi*(59+119*phi))))
          Mr(2,2)=L**2*(16+phi*(123+2*phi*(281+56*phi*(9+10*phi))))
          Mr(3,2)=-(L*(13+phi*(271-16*phi*(4+91*phi))))
          Mr(4,2)=L**2*(5+phi*(45+phi*(-121+56*phi*(-3+10*phi))))
          Mr(5,2)=16*L*(1+phi)**2*(-1+28*phi)
          Mr(1,3)=-4*(-1+8*phi*(47+118*phi))
          Mr(2,3)=-(L*(13+phi*(271-16*phi*(4+91*phi))))
          Mr(3,3)=508+32*phi*(79+134*phi)
          Mr(4,3)=L*(-29+phi*(145+16*phi*(59+119*phi)))
          Mr(5,3)=-512*(1+phi)**2
          Mr(1,4)=-(L*(-13+phi*(-271+16*phi*(4+91*phi))))
          Mr(2,4)=L**2*(5+phi*(45+phi*(-121+56*phi*(-3+10*phi))))
          Mr(3,4)=L*(-29+phi*(145+16*phi*(59+119*phi)))
          Mr(4,4)=L**2*(16+phi*(123+2*phi*(281+56*phi*(9+10*phi))))
          Mr(5,4)=-16*L*(1+phi)**2*(-1+28*phi)
          Mr(1,5)=-512*(1+phi)**2
          Mr(2,5)=16*L*(1+phi)**2*(-1+28*phi)
          Mr(3,5)=-512*(1+phi)**2
          Mr(4,5)=-16*L*(1+phi)**2*(-1+28*phi)
          Mr(5,5)=1024*(1+phi)**2
          Mr=Mr*rho*I/(1+phi)**2/(1+4*phi)**2/210/L
        else
          Mt(1,1)=12*(1046+phi*(14759+2*phi*(40477+660*phi*(157+160*phi))))
          Mt(2,1)=3*L*(228+phi*(2061+5*phi*(1551+4*phi*(997-352*phi*(-3+5*phi)))))
          Mt(3,1)=12*(131+phi*(4139+phi*(30839+2640*phi*(31+25*phi))))
          Mt(4,1)=-3*L*(58+phi*(1591+5*phi*(1463+4*phi*(-7+352*phi*(-3+5*phi)))))
          Mt(5,1)=1056*(1+5*phi)**2*(5+45*phi+84*phi**2)
          Mt(6,1)=-240*L*(1+4*phi)**2*(4+phi*(59-11*phi*(2+5*phi)))
          Mt(1,2)=3*L*(228+phi*(2061+5*phi*(1551+4*phi*(997-352*phi*(-3+5*phi)))))
          Mt(2,2)=L**2*(48+phi*(258+phi*(1129+64*phi*(42+5*phi*(-41+88*phi*(1+10*phi))))))
          Mt(3,2)=3*L*(58+phi*(1591+5*phi*(1463+4*phi*(-7+352*phi*(-3+5*phi)))))
          Mt(4,2)=L**2*(-18+phi*(-402+phi*(-521+64*phi*(42+5*phi*(-41+88*phi*(1+10*phi))))))
          Mt(5,2)=528*L*(1+5*phi)**2*(1+6*phi)
          Mt(6,2)=-8*(L+4*L*phi)**2*(9+phi*(57+10*phi*(-59+44*phi*(11+10*phi))))
          Mt(1,3)=12*(131+phi*(4139+phi*(30839+2640*phi*(31+25*phi))))
          Mt(2,3)=3*L*(58+phi*(1591+5*phi*(1463+4*phi*(-7+352*phi*(-3+5*phi)))))
          Mt(3,3)=12*(1046+phi*(14759+2*phi*(40477+660*phi*(157+160*phi))))
          Mt(4,3)=3*L*(-228+phi*(-2061+5*phi*(-1551+4*phi*(-997+352*phi*(-3+5*phi)))))
          Mt(5,3)=1056*(1+5*phi)**2*(5+45*phi+84*phi**2)
          Mt(6,3)=-240*L*(1+4*phi)**2*(-4+phi*(-59+11*phi*(2+5*phi)))
          Mt(1,4)=-3*L*(58+phi*(1591+5*phi*(1463+4*phi*(-7+352*phi*(-3+5*phi)))))
          Mt(2,4)=L**2*(-18+phi*(-402+phi*(-521+64*phi*(42+5*phi*(-41+88*phi*(1+10*phi))))))
          Mt(3,4)=3*L*(-228+phi*(-2061+5*phi*(-1551+4*phi*(-997+352*phi*(-3+5*phi)))))
          Mt(4,4)=L**2*(48+phi*(258+phi*(1129+64*phi*(42+5*phi*(-41+88*phi*(1+10*phi))))))
          Mt(5,4)=-528*L*(1+5*phi)**2*(1+6*phi)
          Mt(6,4)=-8*(L+4*L*phi)**2*(9+phi*(57+10*phi*(-59+44*phi*(11+10*phi))))
          Mt(1,5)=1056*(1+5*phi)**2*(5+45*phi+84*phi**2)
          Mt(2,5)=528*L*(1+5*phi)**2*(1+6*phi)
          Mt(3,5)=1056*(1+5*phi)**2*(5+45*phi+84*phi**2)
          Mt(4,5)=-528*L*(1+5*phi)**2*(1+6*phi)
          Mt(5,5)=33792*(1+5*phi)**2*(1+3*phi*(3+7*phi))
          Mt(6,5)=0
          Mt(1,6)=-240*L*(1+4*phi)**2*(4+phi*(59-11*phi*(2+5*phi)))
          Mt(2,6)=-8*(L+4*L*phi)**2*(9+phi*(57+10*phi*(-59+44*phi*(11+10*phi))))
          Mt(3,6)=-240*L*(1+4*phi)**2*(-4+phi*(-59+11*phi*(2+5*phi)))
          Mt(4,6)=-8*(L+4*L*phi)**2*(9+phi*(57+10*phi*(-59+44*phi*(11+10*phi))))
          Mt(5,6)=0
          Mt(6,6)=256*(L+4*L*phi)**2*(3+phi*(72+5*phi*(104+11*phi*(13+5*phi))))
          Mt=Mt*rho*A*L/(1+5*phi)**2/(1+4*phi)**2/83160
          Mr(1,1)=12*(139+40*phi*(31+70*phi))
          Mr(2,1)=-3*L*(-13+3*phi*(143+80*phi*(19+45*phi)))
          Mr(3,1)=12*(-11+40*phi*(1+10*phi))
          Mr(4,1)=3*L*(-3+phi*(-141+80*phi*(-6+5*phi)))
          Mr(5,1)=-1536*(1+5*phi)**2
          Mr(6,1)=-240*L*(1+4*phi)**2*(-1+5*phi)
          Mr(1,2)=-3*L*(-13+3*phi*(143+80*phi*(19+45*phi)))
          Mr(2,2)=L**2*(28+phi*(421+2*phi*(2127+40*phi*(281+530*phi))))
          Mr(3,2)=3*L*(3+phi*(141+80*(6-5*phi)*phi))
          Mr(4,2)=L**2*(-5+phi*(-77+phi*(69+40*phi*(37+10*phi))))
          Mr(5,2)=48*L*(1+5*phi)**2*(-1+28*phi)
          Mr(6,2)=8*(L+4*L*phi)**2*(-1+25*phi*(-2+5*phi))
          Mr(1,3)=12*(-11+40*phi*(1+10*phi))
          Mr(2,3)=3*L*(3+phi*(141+80*(6-5*phi)*phi))
          Mr(3,3)=12*(139+40*phi*(31+70*phi))
          Mr(4,3)=3*L*(-13+3*phi*(143+80*phi*(19+45*phi)))
          Mr(5,3)=-1536*(1+5*phi)**2
          Mr(6,3)=240*L*(1+4*phi)**2*(-1+5*phi)
          Mr(1,4)=3*L*(-3+phi*(-141+80*phi*(-6+5*phi)))
          Mr(2,4)=L**2*(-5+phi*(-77+phi*(69+40*phi*(37+10*phi))))
          Mr(3,4)=3*L*(-13+3*phi*(143+80*phi*(19+45*phi)))
          Mr(4,4)=L**2*(28+phi*(421+2*phi*(2127+40*phi*(281+530*phi))))
          Mr(5,4)=-48*L*(1+5*phi)**2*(-1+28*phi)
          Mr(6,4)=8*(L+4*L*phi)**2*(-1+25*phi*(-2+5*phi))
          Mr(1,5)=-1536*(1+5*phi)**2
          Mr(2,5)=48*L*(1+5*phi)**2*(-1+28*phi)
          Mr(3,5)=-1536*(1+5*phi)**2
          Mr(4,5)=-48*L*(1+5*phi)**2*(-1+28*phi)
          Mr(5,5)=3072*(1+5*phi)**2
          Mr(6,5)=0
          Mr(1,6)=-240*L*(1+4*phi)**2*(-1+5*phi)
          Mr(2,6)=8*(L+4*L*phi)**2*(-1+25*phi*(-2+5*phi))
          Mr(3,6)=240*L*(1+4*phi)**2*(-1+5*phi)
          Mr(4,6)=8*(L+4*L*phi)**2*(-1+25*phi*(-2+5*phi))
          Mr(5,6)=0
          Mr(6,6)=256*(L+4*L*phi)**2*(1+5*phi*(1+5*phi))
          Mr=Mr*rho*I/(1+5*phi)**2/(1+4*phi)**2/630/L
        end if
      case default
        call fbem_error_message(error_unit,0,__FILE__,__LINE__,'invalid type of strbeam')
    end select
    fbem_fem_strbeam_Ml=Mt+Mr
  end function fbem_fem_strbeam_Ml

  !! Lateral equivalent forces matrix
  function fbem_fem_strbeam_Ql(etype,esubtype,L,phi)
    implicit none
    integer           :: etype               !! Type of element
    integer           :: esubtype            !! Subtype (for line3):
    real(kind=real64) :: L                   !! Length
    real(kind=real64) :: phi                 !! Ratio of the beam bending stiffness to the shear stiffness
    real(kind=real64) :: fbem_fem_strbeam_Ql(2*fbem_n_nodes(etype),fbem_n_nodes(etype))
    fbem_fem_strbeam_Ql=0
    select case (etype)
      case (fbem_line2)
        fbem_fem_strbeam_Ql(1,1)=42+40*phi
        fbem_fem_strbeam_Ql(2,1)=L*(6+5*phi)
        fbem_fem_strbeam_Ql(3,1)=18+20*phi
        fbem_fem_strbeam_Ql(4,1)=-(L*(4+5*phi))
        fbem_fem_strbeam_Ql(1,2)=18+20*phi
        fbem_fem_strbeam_Ql(2,2)=L*(4+5*phi)
        fbem_fem_strbeam_Ql(3,2)=42+40*phi
        fbem_fem_strbeam_Ql(4,2)=-(L*(6+5*phi))
        fbem_fem_strbeam_Ql=fbem_fem_strbeam_Ql*L/(1+phi)/120
      case (fbem_line3)
        if (esubtype.eq.0) then
          fbem_fem_strbeam_Ql(1,1)=138+628*phi+448*phi**2
          fbem_fem_strbeam_Ql(2,1)=L*(10+31*phi)
          fbem_fem_strbeam_Ql(3,1)=-2*(15+92*phi+56*phi**2)
          fbem_fem_strbeam_Ql(4,1)=L*(4+25*phi)
          fbem_fem_strbeam_Ql(5,1)=32*(1+phi)*(1+7*phi)
          fbem_fem_strbeam_Ql(1,2)=-2*(15+92*phi+56*phi**2)
          fbem_fem_strbeam_Ql(2,2)=-(L*(4+25*phi))
          fbem_fem_strbeam_Ql(3,2)=138+628*phi+448*phi**2
          fbem_fem_strbeam_Ql(4,2)=-(L*(10+31*phi))
          fbem_fem_strbeam_Ql(5,2)=32*(1+phi)*(1+7*phi)
          fbem_fem_strbeam_Ql(1,3)=8*(1+phi)*(11+28*phi)
          fbem_fem_strbeam_Ql(2,3)=8*L*(1+phi)
          fbem_fem_strbeam_Ql(3,3)=8*(1+phi)*(11+28*phi)
          fbem_fem_strbeam_Ql(4,3)=-8*L*(1+phi)
          fbem_fem_strbeam_Ql(5,3)=128*(1+phi)*(3+14*phi)
          fbem_fem_strbeam_Ql=fbem_fem_strbeam_Ql*L/(1+phi)/(1+4*phi)/840
        else
          fbem_fem_strbeam_Ql(1,1)=6*(57+374*phi+560*phi**2)
          fbem_fem_strbeam_Ql(2,1)=L*(18+phi*(21+40*phi*(1+28*phi)))
          fbem_fem_strbeam_Ql(3,1)=6*(-3+8*phi*(8+35*phi))
          fbem_fem_strbeam_Ql(4,1)=L*phi*(-69+40*phi*(1+28*phi))
          fbem_fem_strbeam_Ql(5,1)=96*(1+5*phi)*(1+7*phi)
          fbem_fem_strbeam_Ql(6,1)=-16*L*(1+4*phi)*(3+5*phi*(9+7*phi))
          fbem_fem_strbeam_Ql(1,2)=6*(-3+8*phi*(8+35*phi))
          fbem_fem_strbeam_Ql(2,2)=-(L*phi*(-69+40*phi*(1+28*phi)))
          fbem_fem_strbeam_Ql(3,2)=6*(57+374*phi+560*phi**2)
          fbem_fem_strbeam_Ql(4,2)=-(L*(18+phi*(21+40*phi*(1+28*phi))))
          fbem_fem_strbeam_Ql(5,2)=96*(1+5*phi)*(1+7*phi)
          fbem_fem_strbeam_Ql(6,2)=16*L*(1+4*phi)*(3+5*phi*(9+7*phi))
          fbem_fem_strbeam_Ql(1,3)=24*(1+5*phi)*(11+28*phi)
          fbem_fem_strbeam_Ql(2,3)=24*L*(1+5*phi)
          fbem_fem_strbeam_Ql(3,3)=24*(1+5*phi)*(11+28*phi)
          fbem_fem_strbeam_Ql(4,3)=-24*L*(1+5*phi)
          fbem_fem_strbeam_Ql(5,3)=384*(1+5*phi)*(3+14*phi)
          fbem_fem_strbeam_Ql(6,3)=0
          fbem_fem_strbeam_Ql=fbem_fem_strbeam_Ql*L/(1+5*phi)/(1+4*phi)/2520
        end if
      case default
        call fbem_error_message(error_unit,0,__FILE__,__LINE__,'invalid type of strbeam')
    end select
  end function fbem_fem_strbeam_Ql

  !! Stiffness matrix in local coordinates
  function fbem_fem_strbeam_K(rn,etype,esubtype,theory,L,A,I,k,E,nu)
    implicit none
    integer           :: rn                  !! Ambient space
    integer           :: etype               !! Type of element
    integer           :: esubtype            !! Subtype (for line3)
    integer           :: theory              !! 1: Euler-Bernoulli, 2: Timoshenko
    real(kind=real64) :: L                   !! Length
    real(kind=real64) :: A                   !! Length
    real(kind=real64) :: I(3)                !! Moments of inertia:       I(1)=I11, I(2)=I22, I(3)=I33 (2D analysis).
    real(kind=real64) :: k(3)                !! Shear correction factors: k(1)=kt , k(2)=k2 (2D analysis), k(3)=k3 .
    real(kind=real64) :: E                   !! Young' modulus
    real(kind=real64) :: nu                  !! Poisson's ratio
    real(kind=real64) :: phi2, phi3          !! Ratio of the beam bending stiffness to the shear stiffness
    real(kind=real64) :: G                   !! Shear modulus
    real(kind=real64) :: fbem_fem_strbeam_K(3*(rn-1)*fbem_n_nodes(etype),3*(rn-1)*fbem_n_nodes(etype))
    real(kind=real64) :: Kau(fbem_n_nodes(etype),fbem_n_nodes(etype))
    real(kind=real64) :: Kat(fbem_n_nodes(etype),fbem_n_nodes(etype)), Kat2(2,2)
    real(kind=real64) :: Kl2(2*fbem_n_nodes(etype),2*fbem_n_nodes(etype))
    real(kind=real64) :: Kl3(2*fbem_n_nodes(etype),2*fbem_n_nodes(etype))
    real(kind=real64) :: Kl22(4,4)
    real(kind=real64) :: Kl32(4,4)
    ! Initialize
    fbem_fem_strbeam_K=0
    G=E/(2*(1+nu))
    ! Euler-Bernoulli or Timoshenko
    if (theory.eq.1) then
      phi2=0
      phi3=0
    else
      phi2=12*E*I(3)/(L**2*k(2)*G*A)
      phi3=12*E*I(2)/(L**2*k(3)*G*A)
    end if
    select case (etype)
      case (fbem_line2)
        select case (rn)
          case (2)
            Kau=fbem_fem_strbeam_Ka(etype,L,A,E)
            fbem_fem_strbeam_K([1,4],[1,4])=Kau
            Kl2=fbem_fem_strbeam_Kl(etype,esubtype,L,I(3),E,phi2)
            fbem_fem_strbeam_K([2,3,5,6],[2,3,5,6])=Kl2
          case (3)
            Kau=fbem_fem_strbeam_Ka(etype,L,A,E)
            fbem_fem_strbeam_K([1,7],[1,7])=Kau
            Kat=fbem_fem_strbeam_Ka(etype,L,I(1)*k(1),G)
            fbem_fem_strbeam_K([4,10],[4,10])=Kat
            Kl2=fbem_fem_strbeam_Kl(etype,esubtype,L,I(3),E,phi2)
            fbem_fem_strbeam_K([2,6,8,12],[2,6,8,12])=Kl2
            Kl3=fbem_fem_strbeam_Kl(etype,esubtype,L,I(2),E,phi3)
            Kl3([2,4],:)=-Kl3([2,4],:)
            Kl3(:,[2,4])=-Kl3(:,[2,4])
            fbem_fem_strbeam_K([3,5,9,11],[3,5,9,11])=Kl3
          case default
            call fbem_error_message(error_unit,0,__FILE__,__LINE__,'invalid value of rn')
        end select
      case (fbem_line3)
        select case (rn)
          case (2)
            Kau=fbem_fem_strbeam_Ka(etype,L,A,E)
            fbem_fem_strbeam_K([1,4,7],[1,4,7])=Kau
            Kl2=fbem_fem_strbeam_Kl(etype,esubtype,L,I(3),E,phi2)
            fbem_fem_strbeam_K([2,3,5,6,8,9],[2,3,5,6,8,9])=Kl2
          case (3)
            Kau=fbem_fem_strbeam_Ka(etype,L,A,E)
            fbem_fem_strbeam_K([1,7,13],[1,7,13])=Kau
            if (esubtype.eq.0) then
              Kat2=fbem_fem_strbeam_Ka(fbem_line2,L,I(1)*k(1),G)
              fbem_fem_strbeam_K([4,10],[4,10])=Kat2
            else
              Kat=fbem_fem_strbeam_Ka(etype,L,I(1)*k(1),G)
              fbem_fem_strbeam_K([4,10,16],[4,10,16])=Kat
            end if
            Kl2=fbem_fem_strbeam_Kl(etype,esubtype,L,I(3),E,phi2)
            fbem_fem_strbeam_K([2,6,8,12,14,18],[2,6,8,12,14,18])=Kl2
            Kl3=fbem_fem_strbeam_Kl(etype,esubtype,L,I(2),E,phi3)
            Kl3([2,4,6],:)=-Kl3([2,4,6],:)
            Kl3(:,[2,4,6])=-Kl3(:,[2,4,6])
            fbem_fem_strbeam_K([3,5,9,11,15,17],[3,5,9,11,15,17])=Kl3
          case default
            call fbem_error_message(error_unit,0,__FILE__,__LINE__,'invalid value of rn')
        end select
      case default
        call fbem_error_message(error_unit,0,__FILE__,__LINE__,'invalid type of strbeam')
    end select
  end function fbem_fem_strbeam_K

  !! Mass matrix in local coordinates
  function fbem_fem_strbeam_M(rn,etype,esubtype,theory,L,A,I,k,E,nu,rho)
    implicit none
    integer           :: rn                  !! Ambient space
    integer           :: etype               !! Type of element
    integer           :: esubtype            !! Subtype (for line3)
    integer           :: theory              !! 1: Euler-Bernoulli, 2: Timoshenko
    real(kind=real64) :: L                   !! Length
    real(kind=real64) :: A                   !! Length
    real(kind=real64) :: I(3)                !! Moments of inertia:       I(1)=I11, I(2)=I22, I(3)=I33 (2D analysis).
    real(kind=real64) :: k(3)                !! Shear correction factors: k(1)=kt , k(2)=k2 (2D analysis), k(3)=k3 .
    real(kind=real64) :: E                   !! Young' modulus
    real(kind=real64) :: nu                  !! Poisson's ratio
    real(kind=real64) :: rho                 !! Density
    real(kind=real64) :: phi2, phi3          !! Ratio of the beam bending stiffness to the shear stiffness
    real(kind=real64) :: G                   !! Shear modulus
    real(kind=real64) :: fbem_fem_strbeam_M(3*(rn-1)*fbem_n_nodes(etype),3*(rn-1)*fbem_n_nodes(etype))
    real(kind=real64) :: Mau(fbem_n_nodes(etype),fbem_n_nodes(etype))
    real(kind=real64) :: Mat(fbem_n_nodes(etype),fbem_n_nodes(etype)), Mat2(2,2)
    real(kind=real64) :: Ml2(2*fbem_n_nodes(etype),2*fbem_n_nodes(etype))
    real(kind=real64) :: Ml3(2*fbem_n_nodes(etype),2*fbem_n_nodes(etype))
    real(kind=real64) :: Ml22(4,4)
    real(kind=real64) :: Ml32(4,4)
    ! Initialize
    fbem_fem_strbeam_M=0
    G=E/(2*(1+nu))
    ! Euler-Bernoulli or Timoshenko
    if (theory.eq.1) then
      phi2=0
      phi3=0
    else
      phi2=12*E*I(3)/(L**2*k(2)*G*A)
      phi3=12*E*I(2)/(L**2*k(3)*G*A)
    end if
    select case (etype)
      case (fbem_line2)
        select case (rn)
          case (2)
            Mau=fbem_fem_strbeam_Ma(etype,L,A,rho)
            fbem_fem_strbeam_M([1,4],[1,4])=Mau
            Ml2=fbem_fem_strbeam_Ml(etype,esubtype,L,A,I(3),rho,phi2)
            fbem_fem_strbeam_M([2,3,5,6],[2,3,5,6])=Ml2
          case (3)
            Mau=fbem_fem_strbeam_Ma(etype,L,A,rho)
            fbem_fem_strbeam_M([1,7],[1,7])=Mau
            Mat=fbem_fem_strbeam_Ma(etype,L,I(1),rho)
            fbem_fem_strbeam_M([4,10],[4,10])=Mat
            Ml2=fbem_fem_strbeam_Ml(etype,esubtype,L,A,I(3),rho,phi2)
            fbem_fem_strbeam_M([2,6,8,12],[2,6,8,12])=Ml2
            Ml3=fbem_fem_strbeam_Ml(etype,esubtype,L,A,I(2),rho,phi3)
            Ml3([2,4],:)=-Ml3([2,4],:)
            Ml3(:,[2,4])=-Ml3(:,[2,4])
            fbem_fem_strbeam_M([3,5,9,11],[3,5,9,11])=Ml3
          case default
            call fbem_error_message(error_unit,0,__FILE__,__LINE__,'invalid value of rn')
        end select
      case (fbem_line3)
        select case (rn)
          case (2)
            Mau=fbem_fem_strbeam_Ma(etype,L,A,rho)
            fbem_fem_strbeam_M([1,4,7],[1,4,7])=Mau
            Ml2=fbem_fem_strbeam_Ml(etype,esubtype,L,A,I(3),rho,phi2)
            fbem_fem_strbeam_M([2,3,5,6,8,9],[2,3,5,6,8,9])=Ml2
          case (3)
            Mau=fbem_fem_strbeam_Ma(etype,L,A,rho)
            fbem_fem_strbeam_M([1,7,13],[1,7,13])=Mau
            if (esubtype.eq.0) then
              Mat2=fbem_fem_strbeam_Ma(fbem_line2,L,I(1),rho)
              fbem_fem_strbeam_M([4,10],[4,10])=Mat2
            else
              Mat=fbem_fem_strbeam_Ma(etype,L,I(1),rho)
              fbem_fem_strbeam_M([4,10,16],[4,10,16])=Mat
            end if
            Ml2=fbem_fem_strbeam_Ml(etype,esubtype,L,A,I(3),rho,phi2)
            fbem_fem_strbeam_M([2,6,8,12,14,18],[2,6,8,12,14,18])=Ml2
            Ml3=fbem_fem_strbeam_Ml(etype,esubtype,L,A,I(2),rho,phi3)
            Ml3([2,4,6],:)=-Ml3([2,4,6],:)
            Ml3(:,[2,4,6])=-Ml3(:,[2,4,6])
            fbem_fem_strbeam_M([3,5,9,11,15,17],[3,5,9,11,15,17])=Ml3
          case default
            call fbem_error_message(error_unit,0,__FILE__,__LINE__,'invalid value of rn')
        end select
      case default
        call fbem_error_message(error_unit,0,__FILE__,__LINE__,'invalid type of strbeam')
    end select
  end function fbem_fem_strbeam_M

  !! Aproximate additional mass matrix for flooded elements and elements under hydrodynamic forces in a fluid in local coordinates
  function fbem_fem_strbeam_Madd(rn,etype,esubtype,L,A,rho)
    implicit none
    ! I/O
    integer           :: rn                  !! Ambient space
    integer           :: etype               !! Type of element
    integer           :: esubtype            !! Subtype (for line3)
    real(kind=real64) :: L                   !! Length
    real(kind=real64) :: A(3)                !! Cross-section area considered for each local translational DOF: x',y',z'
    real(kind=real64) :: rho                 !! Fluid density
    real(kind=real64) :: fbem_fem_strbeam_Madd(3*(rn-1)*fbem_n_nodes(etype),3*(rn-1)*fbem_n_nodes(etype)) !! Mass Matrix
    ! Local
    real(kind=real64) :: Mau(fbem_n_nodes(etype),fbem_n_nodes(etype))
    real(kind=real64) :: Ml2(2*fbem_n_nodes(etype),2*fbem_n_nodes(etype))
    real(kind=real64) :: Ml3(2*fbem_n_nodes(etype),2*fbem_n_nodes(etype))
    real(kind=real64) :: Ml22(4,4)
    real(kind=real64) :: Ml32(4,4)
    ! Initialize
    fbem_fem_strbeam_Madd=0
    select case (etype)
      case (fbem_line2)
        select case (rn)
          case (2)
            Mau=fbem_fem_strbeam_Ma(etype,L,A(1),rho)
            fbem_fem_strbeam_Madd([1,4],[1,4])=Mau
            Ml2=fbem_fem_strbeam_Ml(etype,esubtype,L,A(2),0.d0,rho,0.d0)
            fbem_fem_strbeam_Madd([2,3,5,6],[2,3,5,6])=Ml2
          case (3)
            Mau=fbem_fem_strbeam_Ma(etype,L,A(1),rho)
            fbem_fem_strbeam_Madd([1,7],[1,7])=Mau
            Ml2=fbem_fem_strbeam_Ml(etype,esubtype,L,A(2),0.d0,rho,0.d0)
            fbem_fem_strbeam_Madd([2,6,8,12],[2,6,8,12])=Ml2
            Ml3=fbem_fem_strbeam_Ml(etype,esubtype,L,A(3),0.d0,rho,0.d0)
            Ml3([2,4],:)=-Ml3([2,4],:)
            Ml3(:,[2,4])=-Ml3(:,[2,4])
            fbem_fem_strbeam_Madd([3,5,9,11],[3,5,9,11])=Ml3
          case default
            call fbem_error_message(error_unit,0,__FILE__,__LINE__,'invalid value of rn')
        end select
      case (fbem_line3)
        select case (rn)
          case (2)
            Mau=fbem_fem_strbeam_Ma(etype,L,A(1),rho)
            fbem_fem_strbeam_Madd([1,4,7],[1,4,7])=Mau
            Ml2=fbem_fem_strbeam_Ml(etype,esubtype,L,A(2),0.d0,rho,0.d0)
            fbem_fem_strbeam_Madd([2,3,5,6,8,9],[2,3,5,6,8,9])=Ml2
          case (3)
            Mau=fbem_fem_strbeam_Ma(etype,L,A(1),rho)
            fbem_fem_strbeam_Madd([1,7,13],[1,7,13])=Mau
            Ml2=fbem_fem_strbeam_Ml(etype,esubtype,L,A(2),0.d0,rho,0.d0)
            fbem_fem_strbeam_Madd([2,6,8,12,14,18],[2,6,8,12,14,18])=Ml2
            Ml3=fbem_fem_strbeam_Ml(etype,esubtype,L,A(3),0.d0,rho,0.d0)
            Ml3([2,4,6],:)=-Ml3([2,4,6],:)
            Ml3(:,[2,4,6])=-Ml3(:,[2,4,6])
            fbem_fem_strbeam_Madd([3,5,9,11,15,17],[3,5,9,11,15,17])=Ml3
          case default
            call fbem_error_message(error_unit,0,__FILE__,__LINE__,'invalid value of rn')
        end select
      case default
        call fbem_error_message(error_unit,0,__FILE__,__LINE__,'invalid type of strbeam')
    end select
  end function fbem_fem_strbeam_Madd

  !! Equivalent load of distributed mid-line forces in local coordinates
  function fbem_fem_strbeam_Q(rn,etype,esubtype,theory,L,A,I,k,E,nu)
    implicit none
    integer           :: rn                  !! Ambient space
    integer           :: etype               !! Type of element
    integer           :: esubtype            !! Subtype (for line3)
    integer           :: theory              !! 1: Euler-Bernoulli, 2: Timoshenko
    real(kind=real64) :: L                   !! Length
    real(kind=real64) :: A                   !! Length
    real(kind=real64) :: I(3)                !! Moments of inertia:       I(1)=I11, I(2)=I22, I(3)=I33 (2D analysis).
    real(kind=real64) :: k(3)                !! Shear correction factors: k(1)=kt , k(2)=k2 (2D analysis), k(3)=k3 .
    real(kind=real64) :: E                   !! Young' modulus
    real(kind=real64) :: nu                  !! Poisson's ratio
    real(kind=real64) :: phi2, phi3          !! Ratio of the beam bending stiffness to the shear stiffness
    real(kind=real64) :: G                   !! Shear modulus
    real(kind=real64) :: fbem_fem_strbeam_Q(3*(rn-1)*fbem_n_nodes(etype),rn*fbem_n_nodes(etype))
    real(kind=real64) :: Qau(fbem_n_nodes(etype),fbem_n_nodes(etype))
    real(kind=real64) :: Ql2(2*fbem_n_nodes(etype),fbem_n_nodes(etype))
    real(kind=real64) :: Ql3(2*fbem_n_nodes(etype),fbem_n_nodes(etype))
    ! Initialize
    fbem_fem_strbeam_Q=0
    G=E/(2*(1+nu))
    ! Euler-Bernoulli or Timoshenko
    if (theory.eq.1) then
      phi2=0
      phi3=0
    else
      phi2=12*E*I(3)/(L**2*k(2)*G*A)
      phi3=12*E*I(2)/(L**2*k(3)*G*A)
    end if
    select case (etype)
      case (fbem_line2)
        select case (rn)
          case (2)
            Qau=fbem_fem_strbeam_Qa(etype,L)
            fbem_fem_strbeam_Q([1,4],[1,3])=Qau
            Ql2=fbem_fem_strbeam_Ql(etype,esubtype,L,phi2)
            fbem_fem_strbeam_Q([2,3,5,6],[2,4])=Ql2
          case (3)
            Qau=fbem_fem_strbeam_Qa(etype,L)
            fbem_fem_strbeam_Q([1,7],[1,4])=Qau
            Ql2=fbem_fem_strbeam_Ql(etype,esubtype,L,phi2)
            fbem_fem_strbeam_Q([2,6,8,12],[2,5])=Ql2
            Ql3=fbem_fem_strbeam_Ql(etype,esubtype,L,phi3)
            Ql3([2,4],:)=-Ql3([2,4],:)
            fbem_fem_strbeam_Q([3,5,9,11],[3,6])=Ql3
          case default
            call fbem_error_message(error_unit,0,__FILE__,__LINE__,'invalid value of rn')
        end select
      case (fbem_line3)
        select case (rn)
          case (2)
            Qau=fbem_fem_strbeam_Qa(etype,L)
            fbem_fem_strbeam_Q([1,4,7],[1,3,5])=Qau
            Ql2=fbem_fem_strbeam_Ql(etype,esubtype,L,phi2)
            fbem_fem_strbeam_Q([2,3,5,6,8,9],[2,4,6])=Ql2
          case (3)
            Qau=fbem_fem_strbeam_Qa(etype,L)
            fbem_fem_strbeam_Q([1,7,13],[1,4,7])=Qau
            Ql2=fbem_fem_strbeam_Ql(etype,esubtype,L,phi2)
            fbem_fem_strbeam_Q([2,6,8,12,14,18],[2,5,8])=Ql2
            Ql3=fbem_fem_strbeam_Ql(etype,esubtype,L,phi3)
            Ql3([2,4,6],:)=-Ql3([2,4,6],:)
            fbem_fem_strbeam_Q([3,5,9,11,15,17],[3,6,9])=Ql3
          case default
            call fbem_error_message(error_unit,0,__FILE__,__LINE__,'invalid value of rn')
        end select
      case default
        call fbem_error_message(error_unit,0,__FILE__,__LINE__,'invalid type of strbeam')
    end select
  end function fbem_fem_strbeam_Q

  !! Shape function matrix (return displacements in global coordinates)
  function fbem_fem_strbeam_N(rn,etype,esubtype,theory,x,local_axis,A,I,k,E,nu,nodal_axes,xi) result (N)
    implicit none
    ! I/O
    integer           :: rn                                    !! Ambient space
    integer           :: etype                                 !! Type of element
    integer           :: esubtype                              !! Subtype (for line3)
    integer           :: theory                                !! 1: Euler-Bernoulli, 2: Timoshenko
    real(kind=real64) :: x(rn,fbem_n_nodes(etype))             !! Coordinates of the nodes
    real(kind=real64) :: local_axis(rn,rn)                     !! Beam local axis (assumed correct, no checkings performed): axis 1' (axial) local_axis(:,1), axis 2' (lateral) local_axis(:,2), axis 3' (lateral) local_axis(:,3)
    real(kind=real64) :: A                                     !! Length
    real(kind=real64) :: I(3)                                  !! Moments of inertia:       I(1)=I11, I(2)=I22, I(3)=I33 (2D analysis).
    real(kind=real64) :: k(3)                                  !! Shear correction factors: k(1)=kt , k(2)=k2 (2D analysis), k(3)=k3 .
    real(kind=real64) :: E                                     !! Young' modulus
    real(kind=real64) :: nu                                    !! Poisson's ratio
    real(kind=real64) :: nodal_axes(rn,rn,fbem_n_nodes(etype)) !! Axes for each node DOFs nodal_axes(component,axis,node)
    real(kind=real64) :: xi                                    !! Local coordinate
    real(kind=real64) :: N(rn,3*(rn-1)*fbem_n_nodes(etype))
    ! Local
    integer           :: ki
    real(kind=real64) :: Le                  !! Length
    real(kind=real64) :: phi2, phi3          !! Ratio of the beam bending stiffness to the shear stiffness
    real(kind=real64) :: G                   !! Shear modulus
    real(kind=real64) :: Na(fbem_n_nodes(etype))
    real(kind=real64) :: Nl(2*fbem_n_nodes(etype))
    real(kind=real64) :: Nl2(4)
    real(kind=real64) :: L(3*(rn-1)*fbem_n_nodes(etype),3*(rn-1)*fbem_n_nodes(etype))
    real(kind=real64) :: ex(rn,rn), Ln(rn,rn)
    !
    ! Initialize
    !
    N=0
    Le=sqrt(dot_product(x(:,2)-x(:,1),x(:,2)-x(:,1)))
    G=E/(2*(1+nu))
    ! Euler-Bernoulli or Timoshenko
    if (theory.eq.1) then
      phi2=0
      phi3=0
    else
      phi2=12*E*I(3)/(Le**2*k(2)*G*A)
      phi3=12*E*I(2)/(Le**2*k(3)*G*A)
    end if
    ! Shape function matrix, all in local coordinates
    select case (etype)
      case (fbem_line2)
        select case (rn)
          case (2)
            Na=fbem_fem_strbeam_Na(etype,xi)
            N(1,[1,4])=Na
            Nl=fbem_fem_strbeam_Nl(etype,esubtype,Le,phi2,xi)
            N(2,[2,3,5,6])=Nl
          case (3)
            Na=fbem_fem_strbeam_Na(etype,xi)
            N(1,[1,7])=Na
            Nl=fbem_fem_strbeam_Nl(etype,esubtype,Le,phi2,xi)
            N(2,[2,6,8,12])=Nl
            Nl=fbem_fem_strbeam_Nl(etype,esubtype,Le,phi3,xi)
            Nl([2,4])=-Nl([2,4])
            N(3,[3,5,9,11])=Nl
          case default
            call fbem_error_message(error_unit,0,__FILE__,__LINE__,'invalid value of rn')
        end select
      case (fbem_line3)
        select case (rn)
          case (2)
            Na=fbem_fem_strbeam_Na(etype,xi)
            N(1,[1,4,7])=Na
            Nl=fbem_fem_strbeam_Nl(etype,esubtype,Le,phi2,xi)
            N(2,[2,3,5,6,8,9])=Nl
          case (3)
            Na=fbem_fem_strbeam_Na(etype,xi)
            N(1,[1,7,13])=Na
            Nl=fbem_fem_strbeam_Nl(etype,esubtype,Le,phi2,xi)
            N(2,[2,6,8,12,14,18])=Nl
            Nl=fbem_fem_strbeam_Nl(etype,esubtype,Le,phi3,xi)
            Nl([2,4,6])=-Nl([2,4,6])
            N(3,[3,5,9,11,15,17])=Nl
          case default
            call fbem_error_message(error_unit,0,__FILE__,__LINE__,'invalid value of rn')
        end select
      case default
        call fbem_error_message(error_unit,0,__FILE__,__LINE__,'invalid type of strbeam')
    end select
    ! Calculation of change of coordinate matrices
    L=fbem_fem_strbeam_L_element(rn,etype,local_axis,nodal_axes)
    ex=0
    do ki=1,rn
      ex(ki,ki)=1
    end do
    call fbem_coordinate_transformation_L(rn,local_axis,ex,Ln)
    ! Perform change of coordinates
    N=matmul(Ln,matmul(N,transpose(L)))
  end function fbem_fem_strbeam_N

  !! Calculation of all matrices needed for static analysis
  subroutine fbem_fem_strbeam_K_static(rn,etype,esubtype,theory,x,local_axis,A,I,ksh,E,nu,nodal_axes,K)
    implicit none
    ! I/O
    integer           :: rn                                                           !! Number of dimensions of the ambient space
    integer           :: etype                                                        !! Type of element
    integer           :: esubtype                                                     !! Subtype (for line3): 0: without mid-node rotation, 1: with mid-node rotation
    integer           :: theory                                                       !! 1: Euler-Bernoulli, 2: Timoshenko
    real(kind=real64) :: x(rn,fbem_n_nodes(etype))                                    !! Coordinates of the nodes
    real(kind=real64) :: local_axis(rn,rn)                                            !! Beam local axis (assumed correct, no checkings performed): axis 1' (axial) local_axis(:,1), axis 2' (lateral) local_axis(:,2), axis 3' (lateral) local_axis(:,3)
    real(kind=real64) :: A                                                            !! Length
    real(kind=real64) :: I(3)                                                         !! Moments of inertia:       I(1)=I11, I(2)=I22, I(3)=I33 (2D analysis).
    real(kind=real64) :: ksh(3)                                                       !! Shear correction factors: k(1)=kt , k(2)=k2 (2D analysis), k(3)=k3 .
    real(kind=real64) :: E                                                            !! Young' modulus
    real(kind=real64) :: nu                                                           !! Poisson's ratio
    real(kind=real64) :: nodal_axes(rn,rn,fbem_n_nodes(etype))                        !! Axes for each node DOFs nodal_axes(component,axis,node)
    real(kind=real64) :: K(3*(rn-1)*fbem_n_nodes(etype),3*(rn-1)*fbem_n_nodes(etype)) !! Stiffness matrix
    ! Local
    real(kind=real64) :: Le
    real(kind=real64) :: L(3*(rn-1)*fbem_n_nodes(etype),3*(rn-1)*fbem_n_nodes(etype))
    real(kind=real64) :: Kl(3*(rn-1)*fbem_n_nodes(etype),3*(rn-1)*fbem_n_nodes(etype))
    integer :: j
    ! Auxiliary calculations
    Le=sqrt(dot_product(x(:,2)-x(:,1),x(:,2)-x(:,1)))
    ! Calculation of matrices
    L=fbem_fem_strbeam_L_element(rn,etype,local_axis,nodal_axes)
    Kl=fbem_fem_strbeam_K(rn,etype,esubtype,theory,Le,A,I,ksh,E,nu)
    K=matmul(L,matmul(Kl,transpose(L)))
  end subroutine fbem_fem_strbeam_K_static

  !! Calculation of all matrices needed for time harmonic analysis
  subroutine fbem_fem_strbeam_K_harmonic(rn,omega,etype,esubtype,theory,x,local_axis,A,I,ksh,E,nu,rho,nodal_axes,K)
    implicit none
    ! I/O
    integer              :: rn                                                           !! Number of dimensions of the ambient space
    real(kind=real64)    :: omega                                                        !! Circular frequency
    integer              :: etype                                                        !! Type of element
    integer              :: esubtype                                                     !! Subtype (for line3): 0: without mid-node rotation, 1: with mid-node rotation
    integer              :: theory                                                       !! 1: Euler-Bernoulli, 2: Timoshenko
    real(kind=real64)    :: x(rn,fbem_n_nodes(etype))                                    !! Coordinates of the nodes
    real(kind=real64)    :: local_axis(rn,rn)                                            !! Beam local axis (assumed correct, no checkings performed): axis 1' (axial) local_axis(:,1), axis 2' (lateral) local_axis(:,2), axis 3' (lateral) local_axis(:,3)
    real(kind=real64)    :: A                                                            !! Cross-section area
    real(kind=real64)    :: I(3)                                                         !! Moments of inertia:       I(1)=I11, I(2)=I22, I(3)=I33 (2D analysis).
    real(kind=real64)    :: ksh(3)                                                       !! Shear correction factors: k(1)=kt , k(2)=k2 (2D analysis), k(3)=k3 .
    complex(kind=real64) :: E                                                            !! Young' modulus
    real(kind=real64)    :: nu                                                           !! Poisson's ratio
    real(kind=real64)    :: rho                                                          !! Density
    real(kind=real64)    :: nodal_axes(rn,rn,fbem_n_nodes(etype))                        !! Axes for each node DOFs nodal_axes(component,axis,node)
    complex(kind=real64) :: K(3*(rn-1)*fbem_n_nodes(etype),3*(rn-1)*fbem_n_nodes(etype)) !! Stiffness matrix
    ! Local
    real(kind=real64) :: Le
    real(kind=real64) :: L(3*(rn-1)*fbem_n_nodes(etype),3*(rn-1)*fbem_n_nodes(etype))
    real(kind=real64) :: Kl(3*(rn-1)*fbem_n_nodes(etype),3*(rn-1)*fbem_n_nodes(etype))
    real(kind=real64) :: Ml(3*(rn-1)*fbem_n_nodes(etype),3*(rn-1)*fbem_n_nodes(etype))
    ! Auxiliary calculations
    Le=sqrt(dot_product(x(:,2)-x(:,1),x(:,2)-x(:,1)))
    ! Calculation of matrices
    L=fbem_fem_strbeam_L_element(rn,etype,local_axis,nodal_axes)
    Kl=fbem_fem_strbeam_K(rn,etype,esubtype,theory,Le,A,I,ksh,1.d0,nu)
    Ml=fbem_fem_strbeam_M(rn,etype,esubtype,theory,Le,A,I,ksh,1.d0,nu,rho)
    K=matmul(L,matmul(E*Kl-omega**2*Ml,transpose(L)))
  end subroutine fbem_fem_strbeam_K_harmonic

  !! Calculation of all matrices needed for static analysis
  subroutine fbem_fem_strbeam_Q_midline(rn,etype,esubtype,theory,x,local_axis,A,I,ksh,E,nu,nodal_axes,Q)
    implicit none
    ! I/O
    integer           :: rn                                                           !! Number of dimensions of the ambient space
    integer           :: etype                                                        !! Type of element
    integer           :: esubtype                                                     !! Subtype (for line3): 0: without mid-node rotation, 1: with mid-node rotation
    integer           :: theory                                                       !! 1: Euler-Bernoulli, 2: Timoshenko
    real(kind=real64) :: x(rn,fbem_n_nodes(etype))                                    !! Coordinates of the nodes
    real(kind=real64) :: local_axis(rn,rn)                                            !! Beam local axis (assumed correct, no checkings performed): axis 1' (axial) local_axis(:,1), axis 2' (lateral) local_axis(:,2), axis 3' (lateral) local_axis(:,3)
    real(kind=real64) :: A                                                            !! Length
    real(kind=real64) :: I(3)                                                         !! Moments of inertia:       I(1)=I11, I(2)=I22, I(3)=I33 (2D analysis).
    real(kind=real64) :: ksh(3)                                                       !! Shear correction factors: k(1)=kt , k(2)=k2 (2D analysis), k(3)=k3 .
    real(kind=real64) :: E                                                            !! Young' modulus
    real(kind=real64) :: nu                                                           !! Poisson's ratio
    real(kind=real64) :: nodal_axes(rn,rn,fbem_n_nodes(etype))                        !! Axes for each node DOFs nodal_axes(component,axis,node)
    real(kind=real64) :: Q(3*(rn-1)*fbem_n_nodes(etype),rn*fbem_n_nodes(etype))       !! Distributed mid-line force line load matrix: f=Q*t
    ! Local
    real(kind=real64) :: Le
    real(kind=real64) :: L(3*(rn-1)*fbem_n_nodes(etype),3*(rn-1)*fbem_n_nodes(etype))
    real(kind=real64) :: Ql(3*(rn-1)*fbem_n_nodes(etype),rn*fbem_n_nodes(etype))
    real(kind=real64) :: Ll(rn*fbem_n_nodes(etype),rn*fbem_n_nodes(etype))
    ! Auxiliary calculations
    Le=sqrt(dot_product(x(:,2)-x(:,1),x(:,2)-x(:,1)))
    ! Calculation of matrices
    L=fbem_fem_strbeam_L_element(rn,etype,local_axis,nodal_axes)
    Ll=fbem_fem_strbeam_L_load(rn,etype,local_axis,nodal_axes)
    Ql=fbem_fem_strbeam_Q(rn,etype,esubtype,theory,Le,A,I,ksh,E,nu)
    Q=matmul(L,matmul(Ql,transpose(Ll)))
  end subroutine fbem_fem_strbeam_Q_midline

  !! Calculation of transformation matrix of t=Ll*t'
  function fbem_fem_strbeam_L_load(rn,etype,local_axis,nodal_axes)
    implicit none
    ! I/O
    integer           :: rn                                      !! Number of dimensions of the ambient space
    integer           :: etype                                   !! Type of element
    real(kind=real64) :: local_axis(rn,rn)                       !! Beam local axis (assumed correct, no checkings performed): axis 1' (axial) local_axis(:,1), axis 2' (lateral) local_axis(:,2), axis 3' (lateral) local_axis(:,3)
    real(kind=real64) :: nodal_axes(rn,rn,fbem_n_nodes(etype))   !! Axes for each node DOFs nodal_axes(component,axis,node)
    real(kind=real64) :: fbem_fem_strbeam_L_load(rn*fbem_n_nodes(etype),rn*fbem_n_nodes(etype))
    ! Local
    real(kind=real64) :: L(3*(rn-1)*fbem_n_nodes(etype),3*(rn-1)*fbem_n_nodes(etype))
    integer           :: kn, kc, trans(rn*fbem_n_nodes(etype))
    L=fbem_fem_strbeam_L_element(rn,etype,local_axis,nodal_axes)
    do kn=1,fbem_n_nodes(etype)
      do kc=1,rn
        trans(rn*(kn-1)+kc)=3*(rn-1)*(kn-1)+kc
      end do
    end do
    fbem_fem_strbeam_L_load=L(trans,trans)
  end function fbem_fem_strbeam_L_load

 !! Axial D*B product for stress resultant calculation (axial force or torsional moment)
  function fbem_fem_strbeam_DBa(etype,L,A,E,xi)
    implicit none
    integer           :: etype               !! Type of element
    real(kind=real64) :: L                   !! Length
    real(kind=real64) :: A                   !! For axial displacement: Cross section area; For torsion: torsional constant.
    real(kind=real64) :: E                   !! For axial displacement: Young' modulus; For torsion: shear modulus
    real(kind=real64) :: xi                  !! Reference coordinate (0<=xi<=1)
    real(kind=real64) :: fbem_fem_strbeam_DBa(fbem_n_nodes(etype))
    fbem_fem_strbeam_DBa=0
    select case (etype)
      case (fbem_line2)
        fbem_fem_strbeam_DBa(1)=-1
        fbem_fem_strbeam_DBa(2)= 1
      case (fbem_line3)
        fbem_fem_strbeam_DBa(1)=-3+4*xi
        fbem_fem_strbeam_DBa(2)=-1+4*xi
        fbem_fem_strbeam_DBa(3)=4-8*xi
      case default
        call fbem_error_message(error_unit,0,__FILE__,__LINE__,'invalid type of strbeam')
    end select
    fbem_fem_strbeam_DBa=E*A/L*fbem_fem_strbeam_DBa
  end function fbem_fem_strbeam_DBa

  !! Lateral D*B product for stress resultant calculation (shear force and bending moment)
  function fbem_fem_strbeam_DBl(etype,esubtype,L,I,E,phi,xi)
    implicit none
    integer           :: etype               !! Type of element
    integer           :: esubtype            !! Subtype (for line3)
    real(kind=real64) :: L                   !! Length
    real(kind=real64) :: I                   !! Moment of inertia
    real(kind=real64) :: E                   !! Young' modulus
    real(kind=real64) :: phi                 !! Ratio of the beam bending stiffness to the shear stiffness
    real(kind=real64) :: xi                  !! Reference coordinate (0<=xi<=1)
    real(kind=real64) :: fbem_fem_strbeam_DBl(2,2*fbem_n_nodes(etype))
    fbem_fem_strbeam_DBl=0
    select case (etype)
      case (fbem_line2)
        fbem_fem_strbeam_DBl(1,1)=-12/L
        fbem_fem_strbeam_DBl(2,1)=-6+12*xi
        fbem_fem_strbeam_DBl(1,2)=-6
        fbem_fem_strbeam_DBl(2,2)=-(L*(4-6*xi+phi))
        fbem_fem_strbeam_DBl(1,3)=12/L
        fbem_fem_strbeam_DBl(2,3)=6-12*xi
        fbem_fem_strbeam_DBl(1,4)=-6
        fbem_fem_strbeam_DBl(2,4)=L*(-2+6*xi+phi)
        fbem_fem_strbeam_DBl=(E*I/(1+phi)/L**2)*fbem_fem_strbeam_DBl
      case (fbem_line3)
        if (esubtype.eq.0) then
          fbem_fem_strbeam_DBl(1,1)=(12*(16*xi*(1+phi)-3*(3+4*phi)))/L
          fbem_fem_strbeam_DBl(2,1)=-2*(11+20*phi+6*xi*(8*xi*(1+phi)-3*(3+4*phi)))
          fbem_fem_strbeam_DBl(1,2)=-6*(5+8*phi-8*xi*(1+phi))
          fbem_fem_strbeam_DBl(2,2)=-(L*(8+24*xi**2*(1+phi)+phi*(21+4*phi)-6*xi*(5+8*phi)))
          fbem_fem_strbeam_DBl(1,3)=(12*(-7-4*phi+16*xi*(1+phi)))/L
          fbem_fem_strbeam_DBl(2,3)=-2*(5-4*phi+6*xi*(-7-4*phi+8*xi*(1+phi)))
          fbem_fem_strbeam_DBl(1,4)=-6*(-3+8*xi*(1+phi))
          fbem_fem_strbeam_DBl(2,4)=L*(2-18*xi+24*xi**2*(1+phi)+phi*(-3+4*phi))
          fbem_fem_strbeam_DBl(1,5)=(192*(1-2*xi)*(1+phi))/L
          fbem_fem_strbeam_DBl(2,5)=32*(1+6*(-1+xi)*xi)*(1+phi)
          fbem_fem_strbeam_DBl=(E*I/(1+4*phi)/(1+phi)/L**2)*fbem_fem_strbeam_DBl
        else
          fbem_fem_strbeam_DBl(1,1)=(-12*(33+140*phi+8*xi*(-17-70*phi+15*xi*(1+4*phi))))/L
          fbem_fem_strbeam_DBl(2,1)=-2*(23+100*phi+6*xi*(-33-140*phi+4*xi*(17+70*phi-10*xi*(1+4*phi))))
          fbem_fem_strbeam_DBl(1,2)=-2*(39+32*(4-5*phi)*phi+120*xi**2*(1+2*phi-8*phi**2)+24*xi*(-6+5*phi*(-3+8*phi)))
          fbem_fem_strbeam_DBl(2,2)=L*(2*(-6+xi*(39+8*xi*(-9+5*xi)))+(-3+4*xi)*(19+20*xi*(-3+2*xi))*phi-20*(1+16*(-1+xi)*xi*(-1+2*xi))*phi**2)
          fbem_fem_strbeam_DBl(1,3)=(12*(17+60*phi+8*xi*(-13-50*phi+15*xi*(1+4*phi))))/L
          fbem_fem_strbeam_DBl(2,3)=-2*(-7-20*phi+6*xi*(17+60*phi+4*xi*(-13-50*phi+10*xi*(1+4*phi))))
          fbem_fem_strbeam_DBl(1,4)=-2*(15+8*phi+8*(3*xi*(-4+5*xi)+15*xi*(-1+2*xi)*phi-20*(1+6*(-1+xi)*xi)*phi**2))
          fbem_fem_strbeam_DBl(2,4)=L*(-2+2*xi*(15+8*xi*(-6+5*xi))+phi+8*xi*(2+5*xi*(-3+4*xi))*phi+20*(1-16*(-1+xi)*xi*(-1+2*xi))*phi**2)
          fbem_fem_strbeam_DBl(1,5)=(192*(1-2*xi)*(1+5*phi))/L
          fbem_fem_strbeam_DBl(2,5)=32*(1+6*(-1+xi)*xi)*(1+5*phi)
          fbem_fem_strbeam_DBl(1,6)=-32*(1+4*phi)*(6+5*phi+30*(-1+xi)*xi*(1+phi))
          fbem_fem_strbeam_DBl(2,6)=16*L*(-1+2*xi)*(1+4*phi)*(1+10*(-1+xi)*xi*(1+phi))
          fbem_fem_strbeam_DBl=(E*I/(1+4*phi)/(1+5*phi)/L**2)*fbem_fem_strbeam_DBl
        end if
      case default
        call fbem_error_message(error_unit,0,__FILE__,__LINE__,'invalid type of strbeam')
    end select
  end function fbem_fem_strbeam_DBl

  !! D*B product for stress resultant calculation of strbeam
  function fbem_fem_strbeam_DB(rn,etype,esubtype,theory,L,A,I,k,E,nu,xi)
    implicit none
    ! I/O
    integer           :: rn                  !! Ambient space
    integer           :: etype               !! Type of element
    integer           :: esubtype            !! Subtype (for line3)
    integer           :: theory              !! 1: Euler-Bernoulli, 2: Timoshenko
    real(kind=real64) :: L                   !! Length
    real(kind=real64) :: A                   !! Length
    real(kind=real64) :: I(3)                !! Moments of inertia:       I(1)=I11, I(2)=I22, I(3)=I33 (2D analysis).
    real(kind=real64) :: k(3)                !! Shear correction factors: k(1)=kt , k(2)=k2 (2D analysis), k(3)=k3 .
    real(kind=real64) :: E                   !! Young' modulus
    real(kind=real64) :: nu                  !! Poisson's ratio
    real(kind=real64) :: xi                  !! Reference coordinate (0<=xi<=1)
    real(kind=real64) :: fbem_fem_strbeam_DB(3*(rn-1),3*(rn-1)*fbem_n_nodes(etype))
    ! Local
    real(kind=real64) :: phi2, phi3          ! Ratio of the beam bending stiffness to the shear stiffness
    real(kind=real64) :: G                   ! Shear modulus
    real(kind=real64) :: DBau(fbem_n_nodes(etype))
    real(kind=real64) :: DBat(fbem_n_nodes(etype)), DBat2(2)
    real(kind=real64) :: DBl2(2,2*fbem_n_nodes(etype))
    real(kind=real64) :: DBl3(2,2*fbem_n_nodes(etype))
    ! Initialize
    fbem_fem_strbeam_DB=0
    G=E/(2*(1+nu))
    ! Euler-Bernoulli or Timoshenko
    if (theory.eq.1) then
      phi2=0
      phi3=0
    else
      phi2=12*E*I(3)/(L**2*k(2)*G*A)
      phi3=12*E*I(2)/(L**2*k(3)*G*A)
    end if
    select case (etype)
      case (fbem_line2)
        select case (rn)
          case (2)
            DBau=fbem_fem_strbeam_DBa(etype,L,A,E,xi)
            fbem_fem_strbeam_DB(1,[1,4])=DBau
            DBl2=fbem_fem_strbeam_DBl(etype,esubtype,L,I(3),E,phi2,xi)
            fbem_fem_strbeam_DB([2,3],[2,3,5,6])=DBl2
          case (3)
            DBau=fbem_fem_strbeam_DBa(etype,L,A,E,xi)
            fbem_fem_strbeam_DB(1,[1,7])=DBau
            DBat=fbem_fem_strbeam_DBa(etype,L,I(1)*k(1),G,xi)
            fbem_fem_strbeam_DB(4,[4,10])=DBat
            DBl2=fbem_fem_strbeam_DBl(etype,esubtype,L,I(3),E,phi2,xi)
            fbem_fem_strbeam_DB([2,6],[2,6,8,12])=DBl2
            DBl3=fbem_fem_strbeam_DBl(etype,esubtype,L,I(2),E,phi3,xi)
            DBl3(2,:)=-DBl3(2,:)
            DBl3(:,[2,4])=-DBl3(:,[2,4])
            fbem_fem_strbeam_DB([3,5],[3,5,9,11])=DBl3
          case default
            call fbem_error_message(error_unit,0,__FILE__,__LINE__,'invalid value of rn')
        end select
      case (fbem_line3)
        select case (rn)
          case (2)
            DBau=fbem_fem_strbeam_DBa(etype,L,A,E,xi)
            fbem_fem_strbeam_DB(1,[1,4,7])=DBau
            DBl2=fbem_fem_strbeam_DBl(etype,esubtype,L,I(3),E,phi2,xi)
            fbem_fem_strbeam_DB([2,3],[2,3,5,6,8,9])=DBl2
          case (3)
            DBau=fbem_fem_strbeam_DBa(etype,L,A,E,xi)
            fbem_fem_strbeam_DB(1,[1,7,13])=DBau
            if (esubtype.eq.0) then
              DBat2=fbem_fem_strbeam_DBa(fbem_line2,L,I(1)*k(1),G,xi)
              fbem_fem_strbeam_DB(4,[4,10])=DBat2
            else
              DBat=fbem_fem_strbeam_DBa(etype,L,I(1)*k(1),G,xi)
              fbem_fem_strbeam_DB(4,[4,10,16])=DBat
            end if
            DBl2=fbem_fem_strbeam_DBl(etype,esubtype,L,I(3),E,phi2,xi)
            fbem_fem_strbeam_DB([2,6],[2,6,8,12,14,18])=DBl2
            DBl3=fbem_fem_strbeam_DBl(etype,esubtype,L,I(2),E,phi3,xi)
            DBl3(2,:)=-DBl3(2,:)
            DBl3(:,[2,4,6])=-DBl3(:,[2,4,6])
            fbem_fem_strbeam_DB([3,5],[3,5,9,11,15,17])=DBl3
          case default
            call fbem_error_message(error_unit,0,__FILE__,__LINE__,'invalid value of rn')
        end select
      case default
        call fbem_error_message(error_unit,0,__FILE__,__LINE__,'invalid type of strbeam')
    end select
  end function fbem_fem_strbeam_DB

  ! Stress resultants extrapolated from Gauss points. Remember Fsigma must be later multiplied by E (Young modulus)
  subroutine fbem_fem_strbeam_stress_resultants(rn,etype,esubtype,theory,x,local_axis,A,I,ksh,nu,nodal_axes,setype,sedelta,Fsigma)
    implicit none
    ! I/O
    integer           :: rn                                    !! Number of dimensions of the ambient space
    integer           :: etype                                 !! Type of element
    integer           :: esubtype                              !! Subtype (for line3): 0: without mid-node rotation, 1: with mid-node rotation
    integer           :: theory                                !! 1: Euler-Bernoulli, 2: Timoshenko
    real(kind=real64) :: x(rn,fbem_n_nodes(etype))             !! Coordinates of the nodes
    real(kind=real64) :: local_axis(rn,rn)                     !! Beam local axis (assumed correct, no checkings performed): axis 1' (axial) local_axis(:,1), axis 2' (lateral) local_axis(:,2), axis 3' (lateral) local_axis(:,3)
    real(kind=real64) :: A                                     !! Length
    real(kind=real64) :: I(3)                                  !! Moments of inertia:       I(1)=I11, I(2)=I22, I(3)=I33 (2D analysis).
    real(kind=real64) :: ksh(3)                                !! Shear correction factors: k(1)=kt , k(2)=k2 (2D analysis), k(3)=k3 .
    real(kind=real64) :: nu                                    !! Poisson's ratio
    real(kind=real64) :: nodal_axes(rn,rn,fbem_n_nodes(etype)) !! Axes for each node DOFs nodal_axes(component,axis,node)
    integer           :: setype                                !! Stress interpolation: type of interpolation
    real(kind=real64) :: sedelta                               !! Stress interpolation: delta of type of interpolation
    real(kind=real64), allocatable :: Fsigma(:,:,:)            !! Stress resultants matrix at interpolation points (without multiplying by E)
    ! Local
    integer           :: ksp
    real(kind=real64) :: Le
    real(kind=real64) :: L(3*(rn-1)*fbem_n_nodes(etype),3*(rn-1)*fbem_n_nodes(etype))
    real(kind=real64) :: xi
    !
    ! Initialization
    !
    ! Calculate length
    Le=sqrt(dot_product(x(:,2)-x(:,1),x(:,2)-x(:,1)))
    ! Coordinate transformation matrix
    L=fbem_fem_strbeam_L_element(rn,etype,local_axis,nodal_axes)
    ! Select the stress interpolation scheme (Oñate, Structural Analysis with the FEM, Volume 1, 2009)
    if (theory.eq.1) then
      select case (etype)
        case (fbem_line2)
          setype=fbem_line2
          sedelta=0.422649731d0
        case (fbem_line3)
          setype=fbem_line3
          sedelta=0.225403331d0
        case default
          call fbem_error_message(error_unit,0,__FILE__,__LINE__,'etype={line2,line3,line4}')
      end select
    else
      select case (etype)
        case (fbem_line2)
          setype=fbem_line1
          sedelta=0.d0
        case (fbem_line3)
          setype=fbem_line2
          sedelta=0.422649731d0
        case default
          call fbem_error_message(error_unit,0,__FILE__,__LINE__,'etype={line2,line3,line4}')
      end select
    end if
    ! Sin interpolacion desde puntos de Gauss
    !setype=etype
    !sedelta=0
    ! Local stress resultants matrix at interpolation points
    allocate(Fsigma(3*(rn-1),3*(rn-1)*fbem_n_nodes(etype),fbem_n_nodes(setype)))
    Fsigma=0
    !
    ! Loops through sampling points
    !
    do ksp=1,fbem_n_nodes(setype)
      !
      ! Sampling point coordinate
      !
#     define node ksp
#     define etype setype
#     define delta sedelta
#     include <xi_1d_at_node.rc>
#     undef node
#     undef etype
#     undef delta
      !
      ! Convert -1<=xi<=1 to 0<=xi<=1
      !
      xi=0.5d0*(xi+1)
      Fsigma(:,:,ksp)=fbem_fem_strbeam_DB(rn,etype,esubtype,theory,Le,A,I,ksh,1.d0,nu,xi)
      Fsigma(:,:,ksp)=matmul(Fsigma(:,:,ksp),transpose(L))
    end do
  end subroutine fbem_fem_strbeam_stress_resultants

!~   ! Tests
!~   subroutine fbem_fem_strbeam_test()
!~     implicit none
!~     integer, parameter   :: rn=2
!~     integer, parameter   :: etype=fbem_line2                                          !! Type of element
!~     integer           :: esubtype                                                     !! Subtype (for line3): 0: without mid-node rotation, 1: with mid-node rotation
!~     integer           :: theory                                                       !! 0: Euler-Bernoulli, 1: Timoshenko
!~     real(kind=real64) :: x(rn,fbem_n_nodes(etype))                                    !! Coordinates of the nodes
!~     real(kind=real64) :: local_axis(rn,rn)                                            !! Beam local axis (assumed correct, no checkings performed): axis 1' (axial) local_axis(:,1), axis 2' (lateral) local_axis(:,2), axis 3' (lateral) local_axis(:,3)
!~     real(kind=real64) :: A                                                            !! Length
!~     real(kind=real64) :: I(3)                                                         !! Moments of inertia:       2D (I(1)==I33), 3D(I(1)=I11, I(2)=I22, I(3)=I33).
!~     real(kind=real64) :: ksh(3)                                                       !! Shear correction factors: 2D (k(1)==k2 ), 3D(k(1)=kt , k(2)=k2 , k(3)=k3 ).
!~     real(kind=real64) :: E                                                            !! Young' modulus
!~     real(kind=real64) :: nu                                                           !! Poisson's ratio
!~     real(kind=real64) :: nodal_axes(rn,rn,fbem_n_nodes(etype))                        !! Axes for each node DOFs nodal_axes(component,axis,node)
!~     real(kind=real64) :: L(3*(rn-1)*fbem_n_nodes(etype),3*(rn-1)*fbem_n_nodes(etype)) !! Coordinate transformation matrix to nodal axes
!~     real(kind=real64) :: K(3*(rn-1)*fbem_n_nodes(etype),3*(rn-1)*fbem_n_nodes(etype)) !! Stiffness matrix
!~     real(kind=real64) :: Q(3*(rn-1)*fbem_n_nodes(etype),rn*fbem_n_nodes(etype))       !! Distributed mid-line force line load matrix: f=Q*t
!~     integer           :: kn
!~     real(kind=real64) :: F(3*(rn-1))
!~     real(kind=real64) :: matA(3*(rn-1),3*(rn-1)), InvmatA(3*(rn-1),3*(rn-1)), detA
!~     ! Static test
!~     x(:,1)=[0,0]
!~     x(:,2)=[10,0]
!~     nodal_axes(:,1,1)=[1,0]
!~     nodal_axes(:,2,1)=[0,1]
!~     nodal_axes(:,1,2)=[1,0]
!~     nodal_axes(:,2,2)=[0,1]
!~     local_axis(:,1)=[1,0]
!~     local_axis(:,2)=[0,1]
!~     esubtype=0
!~     theory=1
!~     A=1
!~     I=1
!~     ksh=5.d0/6.d0
!~     E=1
!~     nu=0.3
!~     call fbem_fem_strbeam_K_static(rn,etype,esubtype,theory,x,local_axis,A,I,ksh,E,nu,nodal_axes,K)
!~     F=[1,1,0]
!~     matA=K(4:6,4:6)
!~     call fbem_invert_3x3_matrix(matA,InvmatA,detA)
!~     write(*,*) matmul(InvmatA,F)
!~   end subroutine fbem_fem_strbeam_test

  ! ================================================================================================================================
  ! DEGENERATED BEAM
  ! ================================================================================================================================

  !! Calculate the position vector
  function fbem_fem_degbeam3d_x(etype,x_mn,v_mn,t_mn,xi3d)
    implicit none
    ! I/O
    real(kind=real64) :: fbem_fem_degbeam3d_x(3)      !! Position vector
    integer           :: etype                          !! Type of element (displacements interpolation): lineX
    real(kind=real64) :: x_mn(3,fbem_n_nodes(etype))    !! Position vectors of the mid-node.
    real(kind=real64) :: v_mn(3,3,fbem_n_nodes(etype))  !! Director vector at each mid-node.
    real(kind=real64) :: t_mn(3,fbem_n_nodes(etype))    !! Thicknesses at each mid-node.
    real(kind=real64) :: xi3d(3)                        !! Local coordinate
    ! Local
    integer           :: n_mn                           ! Number of mid-nodes
    integer           :: k                              ! Counter
    real(kind=real64) :: xi                             ! Local coordinate xi1
    real(kind=real64) :: aux(10)                        ! Auxiliary variable needed for shape_functions module resources
    real(kind=real64) :: phi(fbem_n_nodes(etype))       ! phi shape functions
    real(kind=real64) :: varphi2(fbem_n_nodes(etype))   ! varphi2 shape functions
    real(kind=real64) :: varphi3(fbem_n_nodes(etype))   ! varphi3 shape functions
    real(kind=real64) :: x(3)                           ! Position vector
    ! Number of mid-nodes
    n_mn=fbem_n_nodes(etype)
    ! Shape functions
    xi=xi3d(1)
#   define delta 0.d0
#   include <phi_1d.rc>
#   undef delta
    varphi2=phi*0.5d0*t_mn(2,:)*xi3d(2)
    varphi3=phi*0.5d0*t_mn(3,:)*xi3d(3)
    ! Calculate position vector x
    x=0.d0
    do k=1,n_mn
      x=x+phi(k)*x_mn(:,k)+varphi2(k)*v_mn(:,2,k)+varphi3(k)*v_mn(:,3,k)
    end do
    fbem_fem_degbeam3d_x=x
  end function fbem_fem_degbeam3d_x

  !! Calculate element coordinate transformation matrix (local -> global): u=Le·u' for all DOF (displacements and rotations) of the element.
  subroutine fbem_fem_degbeam3d_L_element(etype,v_mn,Le)
    implicit none
    ! I/O
    integer           :: etype                                           !! Type of element (displacements interpolation): lineX
    real(kind=real64) :: v_mn(3,3,fbem_n_nodes(etype))                   !! Director vector at each mid-node.
    real(kind=real64) :: Le(6*fbem_n_nodes(etype),6*fbem_n_nodes(etype)) !! Element coordinate transformation matrix
    ! Local
    integer           :: k, is, ie
    Le=0.
    do k=1,fbem_n_nodes(etype)
      ! Displacements
      is=6*k-5
      ie=6*k-3
      Le(is:ie,is:ie)=v_mn(:,:,k)
      ! Rotations
      is=6*k-2
      ie=6*k
      Le(is:ie,is:ie)=v_mn(:,:,k)
    end do
  end subroutine fbem_fem_degbeam3d_L_element

  !! Calculate load coordinate transformation matrix (local -> global): u=Le·u' where each node has 3 DOF in local (u')
  subroutine fbem_fem_degbeam3d_L_load(etype,v_mn,Ll)
    implicit none
    ! I/O
    integer           :: etype                                           !! Type of element (displacements interpolation): lineX
    real(kind=real64) :: v_mn(3,3,fbem_n_nodes(etype))                   !! Director vector at each mid-node.
    real(kind=real64) :: Ll(3*fbem_n_nodes(etype),3*fbem_n_nodes(etype)) !! Load coordinate transformation matrix
    ! Local
    integer           :: k, is, ie
    Ll=0
    do k=1,fbem_n_nodes(etype)
      is=3*k-2
      ie=3*k
      Ll(is:ie,is:ie)=v_mn(:,:,k)
    end do
  end subroutine fbem_fem_degbeam3d_L_load

  !! Calculate the stiffness matrix K for statics, with real E and K
  subroutine fbem_fem_degbeam3d_K(etype,x_mn,v_mn,t_mn,Dp,ngpil,ngpsh,ngpth,K)
    implicit none
    ! I/O
    integer           :: etype                                           !! Element type: line2, line3, line4
    real(kind=real64) :: x_mn(3,fbem_n_nodes(etype))                     !! Nodal position vectors
    real(kind=real64) :: v_mn(3,3,fbem_n_nodes(etype))                   !! Nodal director vectors (must be orthonormal)
    real(kind=real64) :: t_mn(3,fbem_n_nodes(etype))                     !! Thicknesses at each node (only needed in direction v_2 and v_3)
    real(kind=real64) :: Dp(3,3)                                         !! Constitutive matrix (D11: axial, D22: bending, D33: bending)
    integer           :: ngpil                                           !! Integration of in-line contributions
    integer           :: ngpsh                                           !! Integration of shear contributions
    integer           :: ngpth                                           !! Integration of the cross-section
    real(kind=real64) :: K(6*fbem_n_nodes(etype),6*fbem_n_nodes(etype))  !! Stiffness matrix
    ! Local
    integer           :: n_mn                                 ! Number of mid-nodes
    integer           :: kmn                                  ! Counter of mid-nodes
    integer           :: contribution                         ! Contribution of each energy part
    integer           :: ngp                                  ! Number of Gauss points
    integer           :: kxi1, kxi2, kxi3                     ! Integration points counters
    real(kind=real64) :: xi1, xi2, xi3, w1, w2, w3            ! Curvilinear coordinates and quadrature weights
    real(kind=real64) :: aux(10)                              ! Auxiliary variable needed for shape_functions module resources
    real(kind=real64) :: phi(fbem_n_nodes(etype))             ! phi shape functions
    real(kind=real64) :: dphidxi1(fbem_n_nodes(etype))        ! phi shape functions derivatives with respect to xi_1
    real(kind=real64) :: dphidxi2(fbem_n_nodes(etype))        ! phi shape functions derivatives with respect to xi_2
    real(kind=real64) :: dphidxi3(fbem_n_nodes(etype))        ! phi shape functions derivatives with respect to xi_3
    real(kind=real64) :: varphi2(fbem_n_nodes(etype))         ! varphi shape functions (associated with v_2)
    real(kind=real64) :: varphi3(fbem_n_nodes(etype))         ! varphi shape functions (associated with v_3)
    real(kind=real64) :: dvarphi2dxi1(fbem_n_nodes(etype))    ! varphi shape functions derivatives with respect to xi_1
    real(kind=real64) :: dvarphi2dxi2(fbem_n_nodes(etype))    ! varphi shape functions derivatives with respect to xi_2
    real(kind=real64) :: dvarphi2dxi3(fbem_n_nodes(etype))    ! varphi shape functions derivatives with respect to xi_3
    real(kind=real64) :: dvarphi3dxi1(fbem_n_nodes(etype))    ! varphi shape functions derivatives with respect to xi_1
    real(kind=real64) :: dvarphi3dxi2(fbem_n_nodes(etype))    ! varphi shape functions derivatives with respect to xi_2
    real(kind=real64) :: dvarphi3dxi3(fbem_n_nodes(etype))    ! varphi shape functions derivatives with respect to xi_3
    real(kind=real64) :: J(3,3), H(3,3), detJ                 ! Jacobian matrix, its inverse and the jacobian determinant
    real(kind=real64) :: jw                                   ! det(J) * weights
    real(kind=real64) :: T1(3), T2(3), N(3)                   ! Derivatives of position with respect to curvilinear coordinates
    real(kind=real64) :: ep1(3), ep2(3), ep3(3)               ! Local ortogonal axes
    real(kind=real64) :: E(3,6)                               ! E matrix
    real(kind=real64) :: G(6,9)                               ! G matrix
    real(kind=real64) :: M(9,6)                               ! M matrix (shape function matrices derivatives)
    real(kind=real64) :: B(6,6,fbem_n_nodes(etype))           ! B matrix
    real(kind=real64) :: Bp(3,6,fbem_n_nodes(etype))          ! B' matrix
    real(kind=real64) :: Li(6,6)                              ! Li nodal DOF rotation matrix
    integer           :: ki, kis, kie, kj, kjs, kje           ! Counters and nodal DOF limits
    !
    ! Initialization
    !
    K=0
    n_mn=fbem_n_nodes(etype)
    !
    ! Integrate each part (in-line and shear contributions) using different number of gaussian points
    !
    do contribution=1,2
      if (ngpil.ne.ngpsh) then
        select case (contribution)
          ! Integrate in-line contribution with ngpip number of gaussian points
          case (1)
            ngp=ngpil
          ! Integrate shear contributions with ngpsh number of gaussian points
          case (2)
            ngp=ngpsh
        end select
      else
        ! Integrate all together with with ngpip=ngpsh number of gaussian points
        ngp=ngpil
        if (contribution.eq.2) exit
      end if
      !
      ! Loops through integration points
      !
      ! Axial coordinate
      do kxi1=1,ngp
        ! xi_1 and w_1
        xi1=gl11_xi(kxi1,ngp)
        w1=gl11_w(kxi1,ngp)
        !
        ! Axial shape functions and their first derivatives with respect to xi_1, xi_2 and xi_3 at (xi_1,xi_2,xi_3)
        !
#       define xi xi1
#       define dphidxi dphidxi1
#       define delta 0.0d0
#       include <phi_and_dphidxi_1d.rc>
#       undef xi
#       undef dphidxi
#       undef delta
        dphidxi2=0.d0
        dphidxi3=0.d0
        ! Cross-section coordinates
        do kxi2=1,ngpth
          ! xi_2 and w_2
          xi2=gl11_xi(kxi2,ngpth)
          w2=gl11_w(kxi2,ngpth)
          do kxi3=1,ngpth
            ! xi_3 and w_3
            xi3=gl11_xi(kxi3,ngpth)
            w3=gl11_w(kxi3,ngpth)
            !
            ! Cross-section shape functions and their first derivatives with respect to xi_1, xi_2 and xi_3 at (xi_1,xi_2,xi_3)
            !
            ! v_2 direction
            varphi2=phi*0.5d0*xi2*t_mn(2,:)
            dvarphi2dxi1=dphidxi1*0.5d0*xi2*t_mn(2,:)
            dvarphi2dxi2=phi*0.5d0*t_mn(2,:)
            dvarphi2dxi3=0.d0
            ! v_3 direction
            varphi3=phi*0.5d0*xi3*t_mn(3,:)
            dvarphi3dxi1=dphidxi1*0.5d0*xi3*t_mn(3,:)
            dvarphi3dxi2=0.d0
            dvarphi3dxi3=phi*0.5d0*t_mn(3,:)
            !
            ! Calculate Jacobian matrix at (xi_1,xi_2,xi_3)
            !
            J=0.d0
            do kmn=1,n_mn
              J(1,:)=J(1,:)+dphidxi1(kmn)*x_mn(:,kmn)+dvarphi2dxi1(kmn)*v_mn(:,2,kmn)+dvarphi3dxi1(kmn)*v_mn(:,3,kmn)
              J(2,:)=J(2,:)+dphidxi2(kmn)*x_mn(:,kmn)+dvarphi2dxi2(kmn)*v_mn(:,2,kmn)+dvarphi3dxi2(kmn)*v_mn(:,3,kmn)
              J(3,:)=J(3,:)+dphidxi3(kmn)*x_mn(:,kmn)+dvarphi2dxi3(kmn)*v_mn(:,2,kmn)+dvarphi3dxi3(kmn)*v_mn(:,3,kmn)
            end do
            ! Calculate inv(J) and det(J)
            call fbem_invert_3x3_matrix(J,H,detJ)
            !
            ! Calculate local orthogonal system of coordinates (ep_1,ep_2,ep_3) at (xi_1,xi_2,xi_3)
            !
            N =J(1,:)
            T1=J(2,:)
            T2=J(3,:)
            ! ep_1 = n
            ep1=N
            ep1=ep1/sqrt(dot_product(ep1,ep1))
            ! ep_3 =  N x T1 / |N x T1|
            ep3=fbem_cross_product(N,T1)
            ep3=ep3/sqrt(dot_product(ep3,ep3))
            ! ep_2 = ep_3 x ep_1
            ep2=fbem_cross_product(ep3,ep1)
            ep2=ep2/sqrt(dot_product(ep2,ep2))
            ! Global (x) to local (x') tensor transformation matrix
            E=0.d0
            E(1,1)=ep1(1)**2
            E(1,2)=ep1(2)**2
            E(1,3)=ep1(3)**2
            E(1,4)=ep1(1)*ep1(2)
            E(1,5)=ep1(2)*ep1(3)
            E(1,6)=ep1(1)*ep1(3)
            E(2,1)=2.d0*ep1(1)*ep2(1)
            E(2,2)=2.d0*ep1(2)*ep2(2)
            E(2,3)=2.d0*ep1(3)*ep2(3)
            E(3,1)=2.d0*ep1(1)*ep3(1)
            E(3,2)=2.d0*ep1(2)*ep3(2)
            E(3,3)=2.d0*ep1(3)*ep3(3)
            E(2,4)=ep1(1)*ep2(2)+ep1(2)*ep2(1)
            E(2,5)=ep1(2)*ep2(3)+ep1(3)*ep2(2)
            E(2,6)=ep1(1)*ep2(3)+ep1(3)*ep2(1)
            E(3,4)=ep1(1)*ep3(2)+ep1(2)*ep3(1)
            E(3,5)=ep1(2)*ep3(3)+ep1(3)*ep3(2)
            E(3,6)=ep1(1)*ep3(3)+ep1(3)*ep3(1)
            if (ngpil.ne.ngpsh) then
              select case (contribution)
                ! Integrate in-line contribution => Set zero the shear contribution
                case (1)
                  E(2:3,:)=0.d0
                ! Integrate shear contribution => Set zero the in-line contribution
                case (2)
                  E(1,:)=0.d0
              end select
            end if
            ! Derivative transformation matrix for curvilinear to global cartesian tensor transformation
            G=0.d0
            G(1,1)=H(1,1)
            G(2,2)=H(2,1)
            G(3,3)=H(3,1)
            G(4,1)=H(2,1)
            G(4,2)=H(1,1)
            G(5,2)=H(3,1)
            G(5,3)=H(2,1)
            G(6,1)=H(3,1)
            G(6,3)=H(1,1)
            G(1,4)=H(1,2)
            G(2,5)=H(2,2)
            G(3,6)=H(3,2)
            G(4,4)=H(2,2)
            G(4,5)=H(1,2)
            G(5,5)=H(3,2)
            G(5,6)=H(2,2)
            G(6,4)=H(3,2)
            G(6,6)=H(1,2)
            G(1,7)=H(1,3)
            G(2,8)=H(2,3)
            G(3,9)=H(3,3)
            G(4,7)=H(2,3)
            G(4,8)=H(1,3)
            G(5,8)=H(3,3)
            G(5,9)=H(2,3)
            G(6,7)=H(3,3)
            G(6,9)=H(1,3)
            ! Li nodal DOF rotation matrix
            Li=0.d0
            Li(1,1)=1
            Li(2,2)=1
            Li(3,3)=1
            ! Build matrix B for all nodes
            do kmn=1,n_mn
              ! Matrix of derivatives of shape functions matrices with respect to xi1, xi2 and xi3
              M=0.d0
              M(  1,1)= dphidxi1(kmn)
              M(  2,2)= dphidxi1(kmn)
              M(  3,3)= dphidxi1(kmn)
              M(1:3,4)= dvarphi2dxi1(kmn)*v_mn(:,3,kmn)-dvarphi3dxi1(kmn)*v_mn(:,2,kmn)
              M(1:3,5)= dvarphi3dxi1(kmn)*v_mn(:,1,kmn)
              M(1:3,6)=-dvarphi2dxi1(kmn)*v_mn(:,1,kmn)
              M(  4,1)= dphidxi2(kmn)
              M(  5,2)= dphidxi2(kmn)
              M(  6,3)= dphidxi2(kmn)
              M(4:6,4)= dvarphi2dxi2(kmn)*v_mn(:,3,kmn)-dvarphi3dxi2(kmn)*v_mn(:,2,kmn)
              M(4:6,5)= dvarphi3dxi2(kmn)*v_mn(:,1,kmn)
              M(4:6,6)=-dvarphi2dxi2(kmn)*v_mn(:,1,kmn)
              M(  7,1)= dphidxi3(kmn)
              M(  8,2)= dphidxi3(kmn)
              M(  9,3)= dphidxi3(kmn)
              M(7:9,4)= dvarphi2dxi3(kmn)*v_mn(:,3,kmn)-dvarphi3dxi3(kmn)*v_mn(:,2,kmn)
              M(7:9,5)= dvarphi3dxi3(kmn)*v_mn(:,1,kmn)
              M(7:9,6)=-dvarphi2dxi3(kmn)*v_mn(:,1,kmn)
              ! B matrix for each node
              B(:,:,kmn)=matmul(G,M)
              ! B' matrix for each node
              Bp(:,:,kmn)=matmul(E,B(:,:,kmn))
              ! B' matrix for each node with rotations in global coordinates
              Li(4:6,4:6)=v_mn(:,:,kmn)
              Bp(:,:,kmn)=matmul(Bp(:,:,kmn),transpose(Li))
            end do
            ! det(J) * weights
            jw=detJ*w1*w2*w3
            !
            ! Build the element stiffness matrix
            !
            do ki=1,n_mn
              kis=(ki-1)*6+1
              kie=kis+5
              do kj=ki,n_mn
                kjs=(kj-1)*6+1
                kje=kjs+5
                K(kis:kie,kjs:kje)=K(kis:kie,kjs:kje)+matmul(transpose(Bp(:,:,ki)),matmul(Dp,Bp(:,:,kj)))*jw
              end do
            end do
          end do ! xi_3
        end do ! xi_2
      end do ! xi_1
    end do ! Contributions
    ! Apply symmetry of the stiffness matrix
    do ki=1,6*n_mn
      do kj=ki+1,6*n_mn
        K(kj,ki)=K(ki,kj)
      end do
    end do
  end subroutine fbem_fem_degbeam3d_K

  !! Calculate the body load matrix
  subroutine fbem_fem_degbeam3d_Q_body(etype,x_mn,v_mn,t_mn,ngpil,ngpth,Q)
    implicit none
    ! I/O
    integer           :: etype                                           !! Element type: line2, line3, line4
    real(kind=real64) :: x_mn(3,fbem_n_nodes(etype))                     !! Nodal position vectors
    real(kind=real64) :: v_mn(3,3,fbem_n_nodes(etype))                   !! Nodal director vectors (must be orthonormal)
    real(kind=real64) :: t_mn(3,fbem_n_nodes(etype))                     !! Thicknesses at each node (only needed in direction v_2 and v_3)
    integer           :: ngpil                                           !! Integration of in-line contributions
    integer           :: ngpth                                           !! Integration of the cross-section
    real(kind=real64) :: Q(6*fbem_n_nodes(etype),6*fbem_n_nodes(etype))  !! Body load matrix
    ! Local
    integer           :: n_mn                                 ! Number of mid-nodes
    integer           :: kmn                                  ! Counter of mid-nodes
    integer           :: kxi1, kxi2, kxi3                     ! Integration points counters
    real(kind=real64) :: xi1, xi2, xi3, w1, w2, w3            ! Curvilinear coordinates and quadrature weights
    real(kind=real64) :: aux(10)                              ! Auxiliary variable needed for shape_functions module resources
    real(kind=real64) :: phi(fbem_n_nodes(etype))             ! phi shape functions
    real(kind=real64) :: dphidxi1(fbem_n_nodes(etype))        ! phi shape functions derivatives with respect to xi_1
    real(kind=real64) :: dphidxi2(fbem_n_nodes(etype))        ! phi shape functions derivatives with respect to xi_2
    real(kind=real64) :: dphidxi3(fbem_n_nodes(etype))        ! phi shape functions derivatives with respect to xi_3
    real(kind=real64) :: varphi2(fbem_n_nodes(etype))         ! varphi shape functions (associated with v_2)
    real(kind=real64) :: varphi3(fbem_n_nodes(etype))         ! varphi shape functions (associated with v_3)
    real(kind=real64) :: dvarphi2dxi1(fbem_n_nodes(etype))    ! varphi shape functions derivatives with respect to xi_1
    real(kind=real64) :: dvarphi2dxi2(fbem_n_nodes(etype))    ! varphi shape functions derivatives with respect to xi_2
    real(kind=real64) :: dvarphi2dxi3(fbem_n_nodes(etype))    ! varphi shape functions derivatives with respect to xi_3
    real(kind=real64) :: dvarphi3dxi1(fbem_n_nodes(etype))    ! varphi shape functions derivatives with respect to xi_1
    real(kind=real64) :: dvarphi3dxi2(fbem_n_nodes(etype))    ! varphi shape functions derivatives with respect to xi_2
    real(kind=real64) :: dvarphi3dxi3(fbem_n_nodes(etype))    ! varphi shape functions derivatives with respect to xi_3
    real(kind=real64) :: J(3,3), H(3,3), detJ                 ! Jacobian matrix, its inverse and the jacobian determinant
    real(kind=real64) :: jw                                   ! det(J) * weights
    real(kind=real64) :: N(3,6,fbem_n_nodes(etype))           ! N matrix (shape function matrices)
    real(kind=real64) :: Li(6,6)                              ! Li nodal DOF rotation matrix
    integer           :: ki, kis, kie, kj, kjs, kje           ! Counters and nodal DOF limits
    !
    ! Initialization
    !
    Q=0.d0
    n_mn=fbem_n_nodes(etype)
    !
    ! Loops through integration points
    !
    ! Axial coordinate
    do kxi1=1,ngpil
      ! xi_1 and w_1
      xi1=gl11_xi(kxi1,ngpil)
      w1=gl11_w(kxi1,ngpil)
      !
      ! Axial shape functions and their first derivatives with respect to xi_1, xi_2 and xi_3 at (xi_1,xi_2,xi_3)
      !
#     define xi xi1
#     define dphidxi dphidxi1
#     define delta 0.0d0
#     include <phi_and_dphidxi_1d.rc>
#     undef xi
#     undef dphidxi
#     undef delta
      dphidxi2=0.d0
      dphidxi3=0.d0
      ! Cross-section coordinates
      do kxi2=1,ngpth
        ! xi_2 and w_2
        xi2=gl11_xi(kxi2,ngpth)
        w2=gl11_w(kxi2,ngpth)
        do kxi3=1,ngpth
          ! xi_3 and w_3
          xi3=gl11_xi(kxi3,ngpth)
          w3=gl11_w(kxi3,ngpth)
          !
          ! Cross-section shape functions and their first derivatives with respect to xi_1, xi_2 and xi_3 at (xi_1,xi_2,xi_3)
          !
          ! v_2 direction
          varphi2=phi*0.5d0*xi2*t_mn(2,:)
          dvarphi2dxi1=dphidxi1*0.5d0*xi2*t_mn(2,:)
          dvarphi2dxi2=phi*0.5d0*t_mn(2,:)
          dvarphi2dxi3=0.d0
          ! v_3 direction
          varphi3=phi*0.5d0*xi3*t_mn(3,:)
          dvarphi3dxi1=dphidxi1*0.5d0*xi3*t_mn(3,:)
          dvarphi3dxi2=0.d0
          dvarphi3dxi3=phi*0.5d0*t_mn(3,:)
          !
          ! Calculate Jacobian matrix at (xi_1,xi_2,xi_3)
          !
          J=0.d0
          do kmn=1,n_mn
            J(1,:)=J(1,:)+dphidxi1(kmn)*x_mn(:,kmn)+dvarphi2dxi1(kmn)*v_mn(:,2,kmn)+dvarphi3dxi1(kmn)*v_mn(:,3,kmn)
            J(2,:)=J(2,:)+dphidxi2(kmn)*x_mn(:,kmn)+dvarphi2dxi2(kmn)*v_mn(:,2,kmn)+dvarphi3dxi2(kmn)*v_mn(:,3,kmn)
            J(3,:)=J(3,:)+dphidxi3(kmn)*x_mn(:,kmn)+dvarphi2dxi3(kmn)*v_mn(:,2,kmn)+dvarphi3dxi3(kmn)*v_mn(:,3,kmn)
          end do
          ! Calculate inv(J) and det(J)
          call fbem_invert_3x3_matrix(J,H,detJ)
          ! Li nodal DOF rotation matrix
          Li=0.d0
          Li(1,1)=1.d0
          Li(2,2)=1.d0
          Li(3,3)=1.d0
          ! Build N matrix
          N=0.d0
          do kmn=1,n_mn
            ! N matrix for node kmn
            N(  1,1,kmn)= phi(kmn)
            N(  2,2,kmn)= phi(kmn)
            N(  3,3,kmn)= phi(kmn)
            N(1:3,4,kmn)= varphi2(kmn)*v_mn(:,3,kmn)-varphi3(kmn)*v_mn(:,2,kmn)
            N(1:3,5,kmn)= varphi3(kmn)*v_mn(:,1,kmn)
            N(1:3,6,kmn)=-varphi2(kmn)*v_mn(:,1,kmn)
            ! N matrix for node kmn with rotations in global coordinates
            Li(4:6,4:6)=v_mn(:,:,kmn)
            N(:,:,kmn)=matmul(N(:,:,kmn),transpose(Li))
          end do
          ! det(J) * weights
          jw=detJ*w1*w2*w3
          !
          ! Build the load matrix
          !
          do ki=1,n_mn
            kis=(ki-1)*6+1
            kie=kis+5
            do kj=ki,n_mn
              kjs=(kj-1)*6+1
              kje=kjs+5
              Q(kis:kie,kjs:kje)=Q(kis:kie,kjs:kje)+matmul(transpose(N(:,:,ki)),N(:,:,kj))*jw
            end do
          end do
        end do ! xi_3
      end do ! xi_2
    end do ! xi_1
    ! Apply symmetry of the load matrix
    do ki=1,6*n_mn
      do kj=ki+1,6*n_mn
        Q(kj,ki)=Q(ki,kj)
      end do
    end do
  end subroutine fbem_fem_degbeam3d_Q_body

  !! Calculate the line load matrix
  subroutine fbem_fem_degbeam3d_Q_midline(etype,x_mn,ngp,Q)
    implicit none
    ! I/O
    integer           :: etype                                           !! Element type: line2, line3, line4
    real(kind=real64) :: x_mn(3,fbem_n_nodes(etype))                     !! Nodal position vectors
    integer           :: ngp                                             !! Number of integration points
    real(kind=real64) :: Q(6*fbem_n_nodes(etype),3*fbem_n_nodes(etype))  !! Line load matrix
    ! Local
    integer           :: n_mn                                 ! Number of mid-nodes
    integer           :: kmn                                  ! Counter of mid-nodes
    integer           :: kxi1, kxi2, kxi3                     ! Integration points counters
    real(kind=real64) :: xi1, xi2, xi3, w1, w2, w3            ! Curvilinear coordinates and quadrature weights
    real(kind=real64) :: aux(10)                              ! Auxiliary variable needed for shape_functions module resources
    real(kind=real64) :: phi(fbem_n_nodes(etype))             ! phi shape functions
    real(kind=real64) :: dphidxi1(fbem_n_nodes(etype))        ! phi shape functions derivatives with respect to xi_1
    real(kind=real64) :: T(3), jg                             ! Tangent vector and jacobian
    real(kind=real64) :: jw                                   ! jacobian * weights
    real(kind=real64) :: N(3,3,fbem_n_nodes(etype))           ! N matrix (shape function matrices)
    integer           :: ki, kis, kie, kj, kjs, kje           ! Counters and nodal DOF limits
    !
    ! Initialization
    !
    Q=0.d0
    n_mn=fbem_n_nodes(etype)
    !
    ! Loops through integration points
    !
    ! Axial coordinate
    do kxi1=1,ngp
      ! xi_1 and w_1
      xi1=gl11_xi(kxi1,ngp)
      w1=gl11_w(kxi1,ngp)
      !
      ! Axial shape functions and their first derivatives with respect to xi_1 at (xi_1,0,0)
      !
#     define xi xi1
#     define dphidxi dphidxi1
#     define delta 0.0d0
#     include <phi_and_dphidxi_1d.rc>
#     undef xi
#     undef dphidxi
#     undef delta
      T=0.d0
      do kmn=1,n_mn
        T=T+dphidxi1(kmn)*x_mn(:,kmn)
      end do
      jg=sqrt(dot_product(T,T))
      ! Build N matrices
      N=0.d0
      do kmn=1,n_mn
        N(1,1,kmn)=phi(kmn)
        N(2,2,kmn)=phi(kmn)
        N(3,3,kmn)=phi(kmn)
      end do
      ! jacobian * weights
      jw=jg*w1
      ! Build the load matrix
      do ki=1,n_mn
        kis=(ki-1)*6+1
        kie=kis+2
        do kj=1,n_mn
          kjs=(kj-1)*3+1
          kje=kjs+2
          Q(kis:kie,kjs:kje)=Q(kis:kie,kjs:kje)+matmul(transpose(N(:,:,ki)),N(:,:,kj))*jw
        end do
      end do
    end do ! xi_1
  end subroutine fbem_fem_degbeam3d_Q_midline

  subroutine fbem_fem_degbeam_K_static(rn,plane,etype,x_mn,v_mn,t_mn,E,nu,kappa,intmode,ngp,K)
    implicit none
    ! I/O
    integer           :: rn                                                           !! Number of dimensions of the ambient space
    integer           :: plane                                                        !! If rn==2, type of plane analysis: plane strain (1) or stress (2)
    integer           :: etype                                                        !! Element type: line2, line3, line4
    real(kind=real64) :: x_mn(rn,fbem_n_nodes(etype))                                 !! Nodal position vectors
    real(kind=real64) :: v_mn(3,3,fbem_n_nodes(etype))                                !! Nodal director vectors (must be orthonormal)
    real(kind=real64) :: t_mn(3,fbem_n_nodes(etype))                                  !! Thicknesses at each node (only needed in direction v_2 and v_3)
    real(kind=real64) :: E                                                            !! Young' modulus
    real(kind=real64) :: nu                                                           !! Poisson's ratio
    real(kind=real64) :: kappa(3)                                                     !! Shear correction factors: -, ky', kz'
    integer           :: intmode                                                      !! Integration mode: 0 (full), (1) reduced, (2) selective, (3) user-defined
    integer           :: ngp(3)                                                       !! Number of integration points for user-defined mode: (1) membrane, (2) shear (3) thickness
    real(kind=real64) :: K(3*(rn-1)*fbem_n_nodes(etype),3*(rn-1)*fbem_n_nodes(etype)) !! Stiffness matrix
    ! Local
    integer           :: kn
    real(kind=real64) :: x_mn_3d(3,fbem_n_nodes(etype))                    !! Nodal position vectors
    real(kind=real64) :: Dp(3,3)                                           !! Constitutive matrix (D11: axial, D22: shear, D33: shear)
    real(kind=real64) :: K_3d(6*fbem_n_nodes(etype),6*fbem_n_nodes(etype)) !! Stiffness matrix
    integer           :: adof(3*fbem_n_nodes(etype))                       !! Active DOF for building the relevant 2D matrix from 3D
    integer           :: ngpil                                             !! Integration of in-line contributions
    integer           :: ngpsh                                             !! Integration of shear contributions
    integer           :: ngpth                                             !! Integration of the cross-section
    !
    ! Preprocessing required for unified 2D/3D treatment
    !
    select case (rn)
      case (2)
        ! Coordinates
        x_mn_3d=0
        do kn=1,fbem_n_nodes(etype)
          x_mn_3d(1:2,kn)=x_mn(:,kn)
        end do
        ! Local constitutive matrix
        Dp=0
        select case (plane)
          case (1)
            Dp(1,1)=E/(1.0d0-nu**2)
            Dp(2,2)=kappa(2)*E/(1.0d0-nu**2)/(2.d0*(1.d0+nu/(1.0d0-nu)))
            Dp(3,3)=kappa(3)*E/(1.0d0-nu**2)/(2.d0*(1.d0+nu/(1.0d0-nu)))
          case (2)
            Dp(1,1)=E
            Dp(2,2)=kappa(2)*E/(2.d0*(1.d0+nu))
            Dp(3,3)=kappa(3)*E/(2.d0*(1.d0+nu))
          case default
            call fbem_error_message(error_unit,0,__FILE__,__LINE__,'invalid type of plane analysis')
        end select
      case (3)
        ! Coordinates
        x_mn_3d=x_mn
        ! Local constitutive matrix
        Dp=0
        Dp(1,1)=E
        Dp(2,2)=kappa(2)*E/(2.d0*(1.d0+nu))
        Dp(3,3)=kappa(3)*E/(2.d0*(1.d0+nu))
      case default
        call fbem_error_message(error_unit,0,__FILE__,__LINE__,'invalid ambient space')
    end select
    !
    ! Integration mode
    !
    select case (intmode)
      ! Full integration
      case (0)
        ngpil=fbem_n_nodes(etype)
        ngpsh=fbem_n_nodes(etype)
        ngpth=2
      ! Reduced integration
      case (1)
        ngpil=fbem_n_nodes(etype)-1
        ngpsh=fbem_n_nodes(etype)-1
        ngpth=2
      ! Selective integration
      case (2)
        ngpil=fbem_n_nodes(etype)
        ngpsh=fbem_n_nodes(etype)-1
        ngpth=2
      ! User-defined integration
      case (3)
        ngpil=ngp(1)
        ngpsh=ngp(2)
        ngpth=ngp(3)
      case default
        call fbem_error_message(error_unit,0,__FILE__,__LINE__,'invalid integration mode')
    end select
    call fbem_fem_degbeam3d_K(etype,x_mn_3d,v_mn,t_mn,Dp,ngpil,ngpsh,ngpth,K_3d)
    if (rn.eq.2) then
      do kn=1,fbem_n_nodes(etype)
        adof(3*(kn-1)+1)=6*(kn-1)+1
        adof(3*(kn-1)+2)=6*(kn-1)+2
        adof(3*(kn-1)+3)=6*(kn-1)+6
      end do
      K=K_3d(adof,adof)
    else
      K=K_3d
    end if
  end subroutine fbem_fem_degbeam_K_static

  subroutine fbem_fem_degbeam_K_harmonic(rn,omega,plane,etype,x_mn,v_mn,t_mn,E,nu,kappa,rho,K_intmode,K_intngp,M_intmode,M_intngp,K)
    implicit none
    ! I/O
    integer              :: rn                                                           !! Number of dimensions of the ambient space
    real(kind=real64)    :: omega                                                        !! Circular frequency
    integer              :: plane                                                        !! If rn==2, type of plane analysis: plane strain (1) or stress (2)
    integer              :: etype                                                        !! Element type: line2, line3, line4
    real(kind=real64)    :: x_mn(rn,fbem_n_nodes(etype))                                 !! Nodal position vectors
    real(kind=real64)    :: v_mn(3,3,fbem_n_nodes(etype))                                !! Nodal director vectors (must be orthonormal)
    real(kind=real64)    :: t_mn(3,fbem_n_nodes(etype))                                  !! Thicknesses at each node (only needed in direction v_2 and v_3)
    complex(kind=real64) :: E                                                            !! Young' modulus
    real(kind=real64)    :: nu                                                           !! Poisson's ratio
    real(kind=real64)    :: kappa(3)                                                     !! Shear correction factors: -, ky', kz'
    real(kind=real64)    :: rho                                                          !! Density
    integer              :: K_intmode                                                    !! Integration mode (K): 0 (full), (1) reduced, (2) selective, (3) user-defined
    integer              :: K_intngp(3)                                                  !! Number of integration points for user-defined mode (K): (1) membrane, (2) shear (3) thickness
    integer              :: M_intmode                                                    !! Integration mode (M): 0 (full), (1) reduced, (2) selective, (3) user-defined
    integer              :: M_intngp(3)                                                   !! Number of integration points for user-defined mode (M): (1) line (2) thickness
    complex(kind=real64) :: K(3*(rn-1)*fbem_n_nodes(etype),3*(rn-1)*fbem_n_nodes(etype)) !! Stiffness matrix
    ! Local
    integer           :: kn
    real(kind=real64) :: x_mn_3d(3,fbem_n_nodes(etype))                    !! Nodal position vectors
    real(kind=real64) :: Dp(3,3)                                           !! Constitutive matrix (D11: axial, D22: bending, D33: bending)
    real(kind=real64) :: K_3d(6*fbem_n_nodes(etype),6*fbem_n_nodes(etype)) !! Stiffness matrix
    real(kind=real64) :: M_3d(6*fbem_n_nodes(etype),6*fbem_n_nodes(etype)) !! Mass matrix
    integer           :: adof(3*fbem_n_nodes(etype))                       !! Active DOF for building the relevant 2D matrix from 3D
    integer           :: ngpil                                             !! Integration of in-line contributions
    integer           :: ngpsh                                             !! Integration of shear contributions
    integer           :: ngpth                                             !! Integration of the cross-section
    !
    ! Preprocessing required for unified 2D/3D treatment
    !
    select case (rn)
      case (2)
        ! Coordinates
        x_mn_3d=0
        do kn=1,fbem_n_nodes(etype)
          x_mn_3d(1:2,kn)=x_mn(:,kn)
        end do
        ! Local constitutive matrix
        Dp=0
        select case (plane)
          case (1)
            Dp(1,1)=1.d0/(1.0d0-nu**2)
            Dp(2,2)=kappa(2)/(1.0d0-nu**2)/(2.d0*(1.d0+nu/(1.0d0-nu)))
            Dp(3,3)=kappa(3)/(1.0d0-nu**2)/(2.d0*(1.d0+nu/(1.0d0-nu)))
          case (2)
            Dp(1,1)=1.d0
            Dp(2,2)=kappa(2)/(2.d0*(1.d0+nu))
            Dp(3,3)=kappa(3)/(2.d0*(1.d0+nu))
          case default
            call fbem_error_message(error_unit,0,__FILE__,__LINE__,'invalid type of plane analysis')
        end select
      case (3)
        ! Coordinates
        x_mn_3d=x_mn
        ! Local constitutive matrix
        Dp=0
        Dp(1,1)=1.d0
        Dp(2,2)=kappa(2)/(2.d0*(1.d0+nu))
        Dp(3,3)=kappa(3)/(2.d0*(1.d0+nu))
      case default
        call fbem_error_message(error_unit,0,__FILE__,__LINE__,'invalid ambient space')
    end select
    !
    ! Integration mode (K)
    !
    select case (K_intmode)
      ! Full integration
      case (0)
        ngpil=fbem_n_nodes(etype)
        ngpsh=fbem_n_nodes(etype)
        ngpth=2
      ! Reduced integration
      case (1)
        ngpil=fbem_n_nodes(etype)-1
        ngpsh=fbem_n_nodes(etype)-1
        ngpth=2
      ! Selective integration
      case (2)
        ngpil=fbem_n_nodes(etype)
        ngpsh=fbem_n_nodes(etype)-1
        ngpth=2
      ! User-defined integration
      case (3)
        ngpil=K_intngp(1)
        ngpsh=K_intngp(2)
        ngpth=K_intngp(3)
      case default
        call fbem_error_message(error_unit,0,__FILE__,__LINE__,'invalid integration mode')
    end select
    call fbem_fem_degbeam3d_K(etype,x_mn_3d,v_mn,t_mn,Dp,ngpil,ngpsh,ngpth,K_3d)
    !
    ! Integration mode (M)
    !
    select case (M_intmode)
      ! Full integration
      case (0)
        ngpil=fbem_n_nodes(etype)
        ngpth=2
      ! Reduced/selective integration
      case (1,2)
        ngpil=fbem_n_nodes(etype)-1
        ngpth=2
      ! User-defined integration
      case (3)
        ngpil=M_intngp(1)
        ngpth=M_intngp(2)
      case default
        call fbem_error_message(error_unit,0,__FILE__,__LINE__,'invalid integration mode')
    end select
    call fbem_fem_degbeam3d_Q_body(etype,x_mn_3d,v_mn,t_mn,ngpil,ngpth,M_3d)
    if (rn.eq.2) then
      do kn=1,fbem_n_nodes(etype)
        adof(3*(kn-1)+1)=6*(kn-1)+1
        adof(3*(kn-1)+2)=6*(kn-1)+2
        adof(3*(kn-1)+3)=6*(kn-1)+6
      end do
      K=E*K_3d(adof,adof)-omega**2*rho*M_3d(adof,adof)
    else
      K=E*K_3d-omega**2*rho*M_3d
    end if
  end subroutine fbem_fem_degbeam_K_harmonic

  subroutine fbem_fem_degbeam_Q_midline(rn,etype,x_mn,intmode,ngp,Q)
    implicit none
    ! I/O
    integer           :: rn                                                     !! Number of dimensions of the ambient space
    integer           :: etype                                                  !! Element type: line2, line3, line4
    real(kind=real64) :: x_mn(rn,fbem_n_nodes(etype))                           !! Nodal position vectors
    integer           :: intmode                                                !! Integration mode: 0 (full), (1) reduced, (2) selective, (3) user-defined
    integer           :: ngp(3)                                                 !! Number of integration points for user-defined
    real(kind=real64) :: Q(3*(rn-1)*fbem_n_nodes(etype),rn*fbem_n_nodes(etype)) !! Line load matrix
    ! Local
    integer           :: kn
    real(kind=real64) :: x_mn_3d(3,fbem_n_nodes(etype))                    !! Nodal position vectors
    real(kind=real64) :: Q_3d(6*fbem_n_nodes(etype),3*fbem_n_nodes(etype)) !! Load matrix
    integer           :: adof(3*fbem_n_nodes(etype))                       !! Active DOF for building the relevant 2D matrix from 3D
    integer           :: adofc(2*fbem_n_nodes(etype))                      !! Active DOF for building the relevant 2D matrix from 3D
    integer           :: ngpil                                             !! Number of integration points used
    !
    ! Preprocessing required for unified 2D/3D treatment
    !
    select case (rn)
      case (2)
        ! Coordinates, director vectors, and thicknesses
        x_mn_3d=0
        do kn=1,fbem_n_nodes(etype)
          x_mn_3d(1:2,kn)=x_mn(:,kn)
        end do
      case (3)
        ! Coordinates, director vectors, and thicknesses
        x_mn_3d=x_mn
      case default
        call fbem_error_message(error_unit,0,__FILE__,__LINE__,'invalid ambient space')
    end select
    !
    ! Integration mode
    !
    select case (intmode)
      ! Full integration
      case (0)
        ngpil=fbem_n_nodes(etype)
      ! Reduced integration
      case (1,2)
        ngpil=fbem_n_nodes(etype)-1
      ! Custom integration
      case (3)
        ngpil=ngp(1)
      case default
        call fbem_error_message(error_unit,0,__FILE__,__LINE__,'invalid integration mode')
    end select
    call fbem_fem_degbeam3d_Q_midline(etype,x_mn_3d,ngpil,Q_3d)
    if (rn.eq.2) then
      do kn=1,fbem_n_nodes(etype)
        adof (3*(kn-1)+1)=6*(kn-1)+1
        adof (3*(kn-1)+2)=6*(kn-1)+2
        adof (3*(kn-1)+3)=6*(kn-1)+6
        adofc(2*(kn-1)+1)=3*(kn-1)+1
        adofc(2*(kn-1)+2)=3*(kn-1)+2
      end do
      Q=Q_3d(adof,adofc)
    else
      Q=Q_3d
    end if
  end subroutine fbem_fem_degbeam_Q_midline

  subroutine fbem_fem_degbeam_Q_body(rn,etype,x_mn,v_mn,t_mn,intmode,ngp,Q)
    implicit none
    ! I/O
    integer           :: rn                                                           !! Number of dimensions of the ambient space
    integer           :: etype                                                        !! Element type: line2, line3, line4
    real(kind=real64) :: x_mn(rn,fbem_n_nodes(etype))                                 !! Nodal position vectors
    real(kind=real64) :: v_mn(3,3,fbem_n_nodes(etype))                                !! Nodal director vectors (must be orthonormal)
    real(kind=real64) :: t_mn(3,fbem_n_nodes(etype))                                  !! Thicknesses at each node (only needed in direction v_2 and v_3)
    integer           :: intmode                                                      !! Integration mode: 0 (full), (1) reduced, (2) selective, (3) user-defined
    integer           :: ngp(3)                                                       !! Number of integration points for user-defined mode: (1) line (2) thickness
    real(kind=real64) :: Q(3*(rn-1)*fbem_n_nodes(etype),3*(rn-1)*fbem_n_nodes(etype)) !! Load matrix matrix
    ! Local
    integer              :: kn
    real(kind=real64)    :: x_mn_3d(3,fbem_n_nodes(etype))                    !! Nodal position vectors
    real(kind=real64)    :: Q_3d(6*fbem_n_nodes(etype),6*fbem_n_nodes(etype)) !! Mass matrix
    integer              :: adof(3*fbem_n_nodes(etype))                       !! Active DOF for building the relevant 2D matrix from 3D
    integer              :: ngpil                                             !! Integration of in-line contributions
    integer              :: ngpsh                                             !! Integration of shear contributions
    integer              :: ngpth                                             !! Integration of the cross-section
    !
    ! Preprocessing required for unified 2D/3D treatment
    !
    select case (rn)
      case (2)
        ! Coordinates
        x_mn_3d=0
        do kn=1,fbem_n_nodes(etype)
          x_mn_3d(1:2,kn)=x_mn(:,kn)
        end do
      case (3)
        ! Coordinates
        x_mn_3d=x_mn
      case default
        call fbem_error_message(error_unit,0,__FILE__,__LINE__,'invalid ambient space')
    end select
    !
    ! Integration mode
    !
    select case (intmode)
      ! Full integration
      case (0)
        ngpil=fbem_n_nodes(etype)
        ngpth=2
      ! Reduced/selective integration
      case (1,2)
        ngpil=fbem_n_nodes(etype)-1
        ngpth=2
      ! User-defined integration
      case (3)
        ngpil=ngp(1)
        ngpth=ngp(2)
      case default
        call fbem_error_message(error_unit,0,__FILE__,__LINE__,'invalid integration mode')
    end select
    call fbem_fem_degbeam3d_Q_body(etype,x_mn_3d,v_mn,t_mn,ngpil,ngpth,Q_3d)
    if (rn.eq.2) then
      do kn=1,fbem_n_nodes(etype)
        adof(3*(kn-1)+1)=6*(kn-1)+1
        adof(3*(kn-1)+2)=6*(kn-1)+2
        adof(3*(kn-1)+3)=6*(kn-1)+6
      end do
      Q=Q_3d(adof,adof)
    else
      Q=Q_3d
    end if
  end subroutine fbem_fem_degbeam_Q_body

  subroutine fbem_fem_degbeam_L_load(rn,etype,v_mn,Ll)
    implicit none
    ! I/O
    integer           :: rn                                                 !! Number of dimensions of the ambient space
    integer           :: etype                                              !! Element type: line2, line3, line4
    real(kind=real64) :: v_mn(3,3,fbem_n_nodes(etype))                      !! Nodal director vectors (must be orthonormal)
    real(kind=real64) :: Ll(rn*fbem_n_nodes(etype),rn*fbem_n_nodes(etype))  !! Load coordinate transformation matrix
    ! Local
    integer           :: kn
    integer           :: adof(rn*fbem_n_nodes(etype))                        !! Active DOF for building the relevant 2D matrix from 3D
    integer           :: ngpil                                              !! Number of integration points used
    real(kind=real64) :: Ll_3d(3*fbem_n_nodes(etype),3*fbem_n_nodes(etype)) !! Load coordinate transformation matrix
    call fbem_fem_degbeam3d_L_load(etype,v_mn,Ll_3d)
    if (rn.eq.2) then
      do kn=1,fbem_n_nodes(etype)
        adof (2*(kn-1)+1)=3*(kn-1)+1
        adof (2*(kn-1)+2)=3*(kn-1)+2
      end do
      Ll=Ll_3d(adof,adof)
    else
      Ll=Ll_3d
    end if
  end subroutine fbem_fem_degbeam_L_load

  ! Stress resultants extrapolated from Gauss points. Remember Fsigma must be later multiplied by E (Young modulus)
  subroutine fbem_fem_degbeam_stress_resultants(rn,plane,etype,x_mn,v_mn,t_mn,nu,kappa,setype,sedelta,Fsigma)
    implicit none
    ! I/O
    integer           :: rn                                                           !! Number of dimensions of the ambient space
    integer           :: plane                                                        !! If rn==2, type of plane analysis: plane strain (1) or stress (2)
    integer           :: etype                                                        !! Element type: line2, line3, line4
    real(kind=real64) :: x_mn(rn,fbem_n_nodes(etype))                                 !! Nodal position vectors
    real(kind=real64) :: v_mn(3,3,fbem_n_nodes(etype))                                !! Nodal director vectors (must be orthonormal)
    real(kind=real64) :: t_mn(3,fbem_n_nodes(etype))                                  !! Thicknesses at each node (only needed in direction v_2 and v_3)
    real(kind=real64) :: nu                                                           !! Poisson's ratio
    real(kind=real64) :: kappa(3)                                                     !! Shear correction factors: -, ky', kz'
    integer           :: setype                             !! Stress interpolation: type of interpolation
    real(kind=real64) :: sedelta                            !! Stress interpolation: delta of type of interpolation
    real(kind=real64), allocatable :: Fsigma(:,:,:)         !! Stress resultants matrix at interpolation points (without multiplying by E)
    ! Local
    integer           :: kn
    real(kind=real64) :: x_mn_3d(3,fbem_n_nodes(etype))                    !! Nodal position vectors
    real(kind=real64) :: Dp(3,3)                                           !! Constitutive matrix (D11: axial, D22: shear, D33: shear)
    real(kind=real64), allocatable :: Fsigma_3d(:,:,:)                     !! Stress resultants matrix at interpolation points (without multiplying by E)
    integer           :: adof(3*fbem_n_nodes(etype))                       !! Active DOF for building the relevant 2D matrix from 3D
    !
    ! Preprocessing required for unified 2D/3D treatment
    !
    select case (rn)
      case (2)
        ! Coordinates
        x_mn_3d=0
        do kn=1,fbem_n_nodes(etype)
          x_mn_3d(1:2,kn)=x_mn(:,kn)
        end do
        ! Local constitutive matrix
        Dp=0
        select case (plane)
          case (1)
            Dp(1,1)=1.d0/(1.0d0-nu**2)
            Dp(2,2)=kappa(2)/(1.0d0-nu**2)/(2.d0*(1.d0+nu/(1.0d0-nu)))
            Dp(3,3)=kappa(3)/(1.0d0-nu**2)/(2.d0*(1.d0+nu/(1.0d0-nu)))
          case (2)
            Dp(1,1)=1.d0
            Dp(2,2)=kappa(2)/(2.d0*(1.d0+nu))
            Dp(3,3)=kappa(3)/(2.d0*(1.d0+nu))
          case default
            call fbem_error_message(error_unit,0,__FILE__,__LINE__,'invalid type of plane analysis')
        end select
      case (3)
        ! Coordinates
        x_mn_3d=x_mn
        ! Local constitutive matrix
        Dp=0
        Dp(1,1)=1.d0
        Dp(2,2)=kappa(2)/(2.d0*(1.d0+nu))
        Dp(3,3)=kappa(3)/(2.d0*(1.d0+nu))
      case default
        call fbem_error_message(error_unit,0,__FILE__,__LINE__,'invalid ambient space')
    end select
    call fbem_fem_degbeam3d_stress_resultants(etype,x_mn_3d,v_mn,t_mn,Dp,setype,sedelta,Fsigma_3d)
    if (rn.eq.2) then
      allocate(Fsigma(3,3*fbem_n_nodes(etype),fbem_n_nodes(setype)))
      Fsigma=0
      do kn=1,fbem_n_nodes(etype)
        adof(3*(kn-1)+1)=6*(kn-1)+1
        adof(3*(kn-1)+2)=6*(kn-1)+2
        adof(3*(kn-1)+3)=6*(kn-1)+6
      end do
      Fsigma=Fsigma_3d([1,2,6],adof,:)
    else
      allocate(Fsigma(6,6*fbem_n_nodes(etype),fbem_n_nodes(setype)))
      Fsigma=Fsigma_3d
    end if
    deallocate(Fsigma_3d)
  end subroutine fbem_fem_degbeam_stress_resultants

  ! Stress resultants extrapolated from Gauss points. Remember Fsigma must be later multiplied by E (Young modulus)
  subroutine fbem_fem_degbeam3d_stress_resultants(etype,x_mn,v_mn,t_mn,Dp,setype,sedelta,Fsigma)
    implicit none
    ! I/O
    integer           :: etype                              !! Type of element (displacements interpolation): tri3, tri6, quad4, quad8, quad9.
    real(kind=real64) :: x_mn(3,fbem_n_nodes(etype))        !! Position vectors of the mid-plane nodes.
    real(kind=real64) :: v_mn(3,3,fbem_n_nodes(etype))      !! Local axes for each mid-node for the rotation degrees of freedom.
    real(kind=real64) :: t_mn(3,fbem_n_nodes(etype))        !! Thickness of each mid-node in each direction (only v_3 makes sense).
    real(kind=real64) :: Dp(3,3)                            !! Constitutive matrix (D11: axial, D22: bending, D33: bending)
    integer           :: setype                             !! Stress interpolation: type of interpolation
    real(kind=real64) :: sedelta                            !! Stress interpolation: delta of type of interpolation
    real(kind=real64), allocatable :: Fsigma(:,:,:)         !! Stress resultants matrix at interpolation points (without multiplying by E)
    ! Local
    integer           :: n_mn                                 ! Number of mid-nodes
    integer           :: kmn                                  ! Counter of mid-nodes
    integer           :: ksp
    real(kind=real64), allocatable :: sephi(:)                ! Stress interpolation: shape function values
    integer           :: ngpth                                ! Number of Gauss-Legendre integration points for thickness coordinates (xi1,xi2)
    integer           :: kxi1, kxi2, kxi3                     ! Integration points counters
    real(kind=real64) :: xi1, xi2, xi3, w1, w2, w3            ! Curvilinear coordinates and quadrature weights
    real(kind=real64) :: yp, zp, r(3)                         ! Coordinates with respect to cog in local coordinates
    real(kind=real64) :: aux(10)                              ! Auxiliary variable needed for shape_functions module resources
    real(kind=real64) :: phi(fbem_n_nodes(etype))             ! phi shape functions
    real(kind=real64) :: dphidxi1(fbem_n_nodes(etype))        ! phi shape functions derivatives with respect to xi_1
    real(kind=real64) :: dphidxi2(fbem_n_nodes(etype))        ! phi shape functions derivatives with respect to xi_2
    real(kind=real64) :: dphidxi3(fbem_n_nodes(etype))        ! phi shape functions derivatives with respect to xi_3
    real(kind=real64) :: varphi2(fbem_n_nodes(etype))         ! varphi shape functions (associated with v_2)
    real(kind=real64) :: varphi3(fbem_n_nodes(etype))         ! varphi shape functions (associated with v_3)
    real(kind=real64) :: dvarphi2dxi1(fbem_n_nodes(etype))    ! varphi shape functions derivatives with respect to xi_1
    real(kind=real64) :: dvarphi2dxi2(fbem_n_nodes(etype))    ! varphi shape functions derivatives with respect to xi_2
    real(kind=real64) :: dvarphi2dxi3(fbem_n_nodes(etype))    ! varphi shape functions derivatives with respect to xi_3
    real(kind=real64) :: dvarphi3dxi1(fbem_n_nodes(etype))    ! varphi shape functions derivatives with respect to xi_1
    real(kind=real64) :: dvarphi3dxi2(fbem_n_nodes(etype))    ! varphi shape functions derivatives with respect to xi_2
    real(kind=real64) :: dvarphi3dxi3(fbem_n_nodes(etype))    ! varphi shape functions derivatives with respect to xi_3
    real(kind=real64) :: J(3,3), H(3,3), detJ                 ! Jacobian matrix, its inverse and the jacobian determinant
    real(kind=real64) :: jw                                   ! det(J) * weights
    real(kind=real64) :: T1(3), T2(3), N(3)                   ! Derivatives of position with respect to curvilinear coordinates
    real(kind=real64) :: ep1(3), ep2(3), ep3(3)               ! Local ortogonal axes
    real(kind=real64) :: E(3,6)                               ! E matrix
    real(kind=real64) :: G(6,9)                               ! G matrix
    real(kind=real64) :: M(9,6)                               ! M matrix (shape function matrices derivatives)
    real(kind=real64) :: B(6,6,fbem_n_nodes(etype))           ! B matrix
    real(kind=real64) :: Bp(3,6,fbem_n_nodes(etype))          ! B' matrix
    real(kind=real64) :: Li(6,6)                              ! Li nodal DOF rotation matrix
    integer           :: ki, kis, kie, kj, kjs, kje           ! Counters and nodal DOF limits
    !
    ! Initialization
    !
    ! Number of integration points along thickness direction
    ngpth=2
    ! Number of mid-nodes
    n_mn=fbem_n_nodes(etype)
    ! Select the stress interpolation scheme
    select case (etype)
      case (fbem_line2)
        setype=fbem_line1
        sedelta=0.d0
      case (fbem_line3)
        setype=fbem_line2
        sedelta=0.422649731d0
      case (fbem_line4)
        setype=fbem_line3
        sedelta=0.225403331d0
      case default
        call fbem_error_message(error_unit,0,__FILE__,__LINE__,'etype={line2,line3,line4}')
    end select
    ! Sin interpolacion desde puntos de Gauss
    !setype=etype
    !sedelta=0
    ! Local stress resultants matrix at interpolation points
    allocate(Fsigma(6,6*n_mn,fbem_n_nodes(setype)))
    Fsigma=0
    !
    ! Loops through sampling points
    !
    do ksp=1,fbem_n_nodes(setype)
      !
      ! Sampling point coordinate
      !
#     define node ksp
#     define etype setype
#     define delta sedelta
#     define xi xi1
#     include <xi_1d_at_node.rc>
#     undef node
#     undef etype
#     undef delta
#     undef xi
      !
      ! Axial shape functions and their first derivatives with respect to xi_1, xi_2 and xi_3 at (xi_1,xi_2,xi_3)
      !
#     define xi xi1
#     define dphidxi dphidxi1
#     define delta 0.0d0
#     include <phi_and_dphidxi_1d.rc>
#     undef xi
#     undef dphidxi
#     undef delta
      dphidxi2=0.d0
      dphidxi3=0.d0
      ! Cross-section coordinates
      do kxi2=1,ngpth
        ! xi_2 and w_2
        xi2=gl11_xi(kxi2,ngpth)
        w2=gl11_w(kxi2,ngpth)
        do kxi3=1,ngpth
          ! xi_3 and w_3
          xi3=gl11_xi(kxi3,ngpth)
          w3=gl11_w(kxi3,ngpth)
          !
          ! Cross-section shape functions and their first derivatives with respect to xi_1, xi_2 and xi_3 at (xi_1,xi_2,xi_3)
          !
          ! v_2 direction
          varphi2=phi*0.5d0*xi2*t_mn(2,:)
          dvarphi2dxi1=dphidxi1*0.5d0*xi2*t_mn(2,:)
          dvarphi2dxi2=phi*0.5d0*t_mn(2,:)
          dvarphi2dxi3=0.d0
          ! v_3 direction
          varphi3=phi*0.5d0*xi3*t_mn(3,:)
          dvarphi3dxi1=dphidxi1*0.5d0*xi3*t_mn(3,:)
          dvarphi3dxi2=0.d0
          dvarphi3dxi3=phi*0.5d0*t_mn(3,:)
          !
          ! Calculate Jacobian matrix at (xi_1,xi_2,xi_3)
          !
          J=0.d0
          do kmn=1,n_mn
            J(1,:)=J(1,:)+dphidxi1(kmn)*x_mn(:,kmn)+dvarphi2dxi1(kmn)*v_mn(:,2,kmn)+dvarphi3dxi1(kmn)*v_mn(:,3,kmn)
            J(2,:)=J(2,:)+dphidxi2(kmn)*x_mn(:,kmn)+dvarphi2dxi2(kmn)*v_mn(:,2,kmn)+dvarphi3dxi2(kmn)*v_mn(:,3,kmn)
            J(3,:)=J(3,:)+dphidxi3(kmn)*x_mn(:,kmn)+dvarphi2dxi3(kmn)*v_mn(:,2,kmn)+dvarphi3dxi3(kmn)*v_mn(:,3,kmn)
          end do
          ! Calculate inv(J) and det(J)
          call fbem_invert_3x3_matrix(J,H,detJ)
          ! Calculate local orthogonal system of coordinates (ep_1,ep_2,ep_3) at (xi_1,xi_2,xi_3)
          N =J(1,:)
          T1=J(2,:)
          T2=J(3,:)
          call fbem_fem_degbeam3d_stress_resultants_ep(N,T1,T2,ep1,ep2,ep3)
          ! Global (x) to local (x') tensor transformation matrix
          E=0.d0
          E(1,1)=ep1(1)**2
          E(1,2)=ep1(2)**2
          E(1,3)=ep1(3)**2
          E(1,4)=ep1(1)*ep1(2)
          E(1,5)=ep1(2)*ep1(3)
          E(1,6)=ep1(1)*ep1(3)
          E(2,1)=2.d0*ep1(1)*ep2(1)
          E(2,2)=2.d0*ep1(2)*ep2(2)
          E(2,3)=2.d0*ep1(3)*ep2(3)
          E(3,1)=2.d0*ep1(1)*ep3(1)
          E(3,2)=2.d0*ep1(2)*ep3(2)
          E(3,3)=2.d0*ep1(3)*ep3(3)
          E(2,4)=ep1(1)*ep2(2)+ep1(2)*ep2(1)
          E(2,5)=ep1(2)*ep2(3)+ep1(3)*ep2(2)
          E(2,6)=ep1(1)*ep2(3)+ep1(3)*ep2(1)
          E(3,4)=ep1(1)*ep3(2)+ep1(2)*ep3(1)
          E(3,5)=ep1(2)*ep3(3)+ep1(3)*ep3(2)
          E(3,6)=ep1(1)*ep3(3)+ep1(3)*ep3(1)
          ! Derivative transformation matrix for curvilinear to global cartesian tensor transformation
          G=0.d0
          G(1,1)=H(1,1)
          G(2,2)=H(2,1)
          G(3,3)=H(3,1)
          G(4,1)=H(2,1)
          G(4,2)=H(1,1)
          G(5,2)=H(3,1)
          G(5,3)=H(2,1)
          G(6,1)=H(3,1)
          G(6,3)=H(1,1)
          G(1,4)=H(1,2)
          G(2,5)=H(2,2)
          G(3,6)=H(3,2)
          G(4,4)=H(2,2)
          G(4,5)=H(1,2)
          G(5,5)=H(3,2)
          G(5,6)=H(2,2)
          G(6,4)=H(3,2)
          G(6,6)=H(1,2)
          G(1,7)=H(1,3)
          G(2,8)=H(2,3)
          G(3,9)=H(3,3)
          G(4,7)=H(2,3)
          G(4,8)=H(1,3)
          G(5,8)=H(3,3)
          G(5,9)=H(2,3)
          G(6,7)=H(3,3)
          G(6,9)=H(1,3)
          ! Li nodal DOF rotation matrix
          Li=0.d0
          Li(1,1)=1.d0
          Li(2,2)=1.d0
          Li(3,3)=1.d0
          ! Build matrix B for all nodes
          do kmn=1,n_mn
            ! Matrix of derivatives of shape functions matrices with respect to xi1, xi2 and xi3
            M=0.d0
            M(  1,1)= dphidxi1(kmn)
            M(  2,2)= dphidxi1(kmn)
            M(  3,3)= dphidxi1(kmn)
            M(1:3,4)= dvarphi2dxi1(kmn)*v_mn(:,3,kmn)-dvarphi3dxi1(kmn)*v_mn(:,2,kmn)
            M(1:3,5)= dvarphi3dxi1(kmn)*v_mn(:,1,kmn)
            M(1:3,6)=-dvarphi2dxi1(kmn)*v_mn(:,1,kmn)
            M(  4,1)= dphidxi2(kmn)
            M(  5,2)= dphidxi2(kmn)
            M(  6,3)= dphidxi2(kmn)
            M(4:6,4)= dvarphi2dxi2(kmn)*v_mn(:,3,kmn)-dvarphi3dxi2(kmn)*v_mn(:,2,kmn)
            M(4:6,5)= dvarphi3dxi2(kmn)*v_mn(:,1,kmn)
            M(4:6,6)=-dvarphi2dxi2(kmn)*v_mn(:,1,kmn)
            M(  7,1)= dphidxi3(kmn)
            M(  8,2)= dphidxi3(kmn)
            M(  9,3)= dphidxi3(kmn)
            M(7:9,4)= dvarphi2dxi3(kmn)*v_mn(:,3,kmn)-dvarphi3dxi3(kmn)*v_mn(:,2,kmn)
            M(7:9,5)= dvarphi3dxi3(kmn)*v_mn(:,1,kmn)
            M(7:9,6)=-dvarphi2dxi3(kmn)*v_mn(:,1,kmn)
            ! B matrix for each node
            B(:,:,kmn)=matmul(G,M)
            ! B' matrix for each node
            Bp(:,:,kmn)=matmul(E,B(:,:,kmn))
            ! B' matrix for each node with rotations in global coordinates
            Li(4:6,4:6)=v_mn(:,:,kmn)
            Bp(:,:,kmn)=matmul(Bp(:,:,kmn),transpose(Li))
          end do
          ! |T1xT2| * weights
          N=fbem_cross_product(T1,T2)
          jw=sqrt(dot_product(N,N))*w2*w3
          ! Distances y' and z'
          ! Distance vector between (xi1,0,0)->(xi1,xi2,xi3)
          r=0
          do kmn=1,n_mn
            r=r+varphi2(kmn)*v_mn(:,2,kmn)+varphi3(kmn)*v_mn(:,3,kmn)
          end do
          ! Projection
          yp=dot_product(r,ep2)
          zp=dot_product(r,ep3)
          ! Build the stress resultants matrix at (xi1)
          do kj=1,n_mn
            kjs=(kj-1)*6+1
            kje=kjs+5
            Fsigma(1:3,kjs:kje,ksp)=Fsigma(1:3,kjs:kje,ksp)+matmul(Dp(1:3,:),Bp(:,:,kj))*jw    ! Axial and shear forces: Nx', Vy', Vz'
            Fsigma(  4,kjs:kje,ksp)=Fsigma(  4,kjs:kje,ksp)-matmul(Dp(  2,:),Bp(:,:,kj))*zp*jw ! Torsional moment: Mx' (sigma_x'y' contribution)
            Fsigma(  4,kjs:kje,ksp)=Fsigma(  4,kjs:kje,ksp)+matmul(Dp(  3,:),Bp(:,:,kj))*yp*jw ! Torsional moment: Mx' (sigma_x'z' contribution)
            Fsigma(  5,kjs:kje,ksp)=Fsigma(  5,kjs:kje,ksp)+matmul(Dp(  1,:),Bp(:,:,kj))*zp*jw ! Bending moment: My'
            Fsigma(  6,kjs:kje,ksp)=Fsigma(  6,kjs:kje,ksp)-matmul(Dp(  1,:),Bp(:,:,kj))*yp*jw ! Bending moment: Mz'
          end do
        end do ! xi_3
      end do ! xi_2
    end do ! xi_1 sampling point
  end subroutine fbem_fem_degbeam3d_stress_resultants

  !! Build the orthonormal basis e1', e2', e3' for stress resultant calculation
  subroutine fbem_fem_degbeam3d_stress_resultants_ep(N,T1,T2,ep1,ep2,ep3)
    implicit none
    real(kind=real64) :: N(3), T1(3), T2(3)
    real(kind=real64) :: ep1(3), ep2(3), ep3(3)
    ! ep_1 = n
    ep1=N
    ep1=ep1/sqrt(dot_product(ep1,ep1))
    ! ep_3 =  N x T1 / |N x T1|
    ep3(1)=N(2)*T1(3)-N(3)*T1(2)
    ep3(2)=N(3)*T1(1)-N(1)*T1(3)
    ep3(3)=N(1)*T1(2)-N(2)*T1(1)
    ep3=ep3/sqrt(dot_product(ep3,ep3))
    ! ep_2 = ep_3 x ep_1
    ep2(1)=ep3(2)*ep1(3)-ep3(3)*ep1(2)
    ep2(2)=ep3(3)*ep1(1)-ep3(1)*ep1(3)
    ep2(3)=ep3(1)*ep1(2)-ep3(2)*ep1(1)
    ep2=ep2/sqrt(dot_product(ep2,ep2))
  end subroutine fbem_fem_degbeam3d_stress_resultants_ep

!  subroutine test_3dbeamdeg()
!    implicit none
!    integer :: i, j, k
!    !
!    ! Ejemplo Oñate Vol: curved cantilever beam
!    !
!    ! Elemento cuadratico (3 nodos)
!    real(kind=real64) :: x_md(3,3)
!    real(kind=real64) :: v_md(3,3,3)
!    real(kind=real64) :: t_md(3,3), Dp(3,3)
!    real(kind=real64) :: xi(3), x(3), Ks(6*3,6*3)
!    !
!    x_md (:,1)=4*[0.,5.,0.]
!    v_md(:,1,1)=[1.,0.,0.]
!    t_md(  1,1)=0.
!    v_md(:,2,1)=[0.,1.,0.]
!    t_md(  2,1)=1.
!    v_md(:,3,1)=[0.,0.,1.]
!    t_md(  3,1)=1.
!    !
!    x_md (:,2)=4*[5.,0.,0.]
!    v_md(:,1,2)=[0.,-1.,0.]
!    t_md(  1,2)=0.
!    v_md(:,2,2)=[1.,0.,0.]
!    t_md(  2,2)=1.
!    v_md(:,3,2)=[0.,0.,1.]
!    t_md(  3,2)=1.
!    !
!    x_md (:,3)=4*[3.535533906,3.535533906,0.]
!    v_md(:,1,3)=[0.7071,-0.7071,0.]
!    t_md(  1,3)=0.
!    v_md(:,2,3)=[0.7071,0.7071,0.]
!    t_md(  2,3)=1.
!    v_md(:,3,3)=[0.,0.,1.]
!    t_md(  3,3)=1.
!    !
!    Dp=0.d0
!    Dp(1,1)=1.e6
!    Dp(2,2)=0.41667e6
!    Dp(3,3)=0.41667e6
!    !
!    do i=0,100
!    do j=0,10
!    do k=0,10
!      xi(1)=-1.d0+2.d0*dble(i)/dble(100)
!      xi(2)=-1.d0+2.d0*dble(j)/dble(10)
!      xi(3)=-1.d0+2.d0*dble(k)/dble(10)
!      x=fbem_fem_degbeam3d_x(fbem_line3,x_md,v_md,t_md,xi)
!      write(63,'(6e25.16)') xi, x
!    end do
!    end do
!    end do
!    !
!    call fbem_fem_degbeam3d_K(fbem_line3,x_md,v_md,t_md,Dp,2,2,2,Ks)
!    do i=1,18
!      write(99,'(18e25.16)') Ks(:,i)
!    end do
!!    ! Ejecutar en octave:
!!    S=load("fort.99");
!!    K=S(7:18,7:18);
!!    b=zeros(12,1);
!!    b(3)=1;
!!    x=linsolve(K,b);
!!    x(3)
!!    ! =========================
!!    ! Elemento cubico (4 nodos)
!!    ! =========================
!!    real(kind=real64) :: x_md(3,4)
!!    real(kind=real64) :: v_md(3,3,4)
!!    real(kind=real64) :: t_md(3,4), Dp(3,3)
!!    real(kind=real64) :: xi(3), x(3), Ks(6*4,6*4), Qs(6*4,6*4)
!!    !
!!    x_md (:,1)=4*[0.,5.,0.]
!!    v_md(:,1,1)=[1.,0.,0.]
!!    t_md(  1,1)=0.
!!    v_md(:,2,1)=[0.,1.,0.]
!!    t_md(  2,1)=1.
!!    v_md(:,3,1)=[0.,0.,1.]
!!    t_md(  3,1)=1.
!!    !
!!    x_md (:,2)=4*[5.,0.,0.]
!!    v_md(:,1,2)=[0.,-1.,0.]
!!    t_md(  1,2)=0.
!!    v_md(:,2,2)=[1.,0.,0.]
!!    t_md(  2,2)=1.
!!    v_md(:,3,2)=[0.,0.,1.]
!!    t_md(  3,2)=1.
!!    !
!!    x_md (:,3)=4*[2.5,4.330127019,0.]
!!    v_md(:,1,3)=[0.866025404,-0.5,0.]
!!    t_md(  1,3)=0.
!!    v_md(:,2,3)=[0.5,0.866025404,0.]
!!    t_md(  2,3)=1.
!!    v_md(:,3,3)=[0.,0.,1.]
!!    t_md(  3,3)=1.
!!    !
!!    x_md (:,4)=4*[4.330127019,2.5,0.]
!!    v_md(:,1,4)=[0.5,-0.866025404,0.]
!!    t_md(  1,4)=0.
!!    v_md(:,2,4)=[0.866025404,0.5,0.]
!!    t_md(  2,4)=1.
!!    v_md(:,3,4)=[0.,0.,1.]
!!    t_md(  3,4)=1.
!!    !
!!    Dp=0.d0
!!    Dp(1,1)=1.e6
!!    Dp(2,2)=0.41667e6
!!    Dp(3,3)=0.41667e6
!!    !
!!    do i=0,100
!!    do j=0,10
!!    do k=0,10
!!      xi(1)=-1.d0+2.d0*dble(i)/dble(100)
!!      xi(2)=-1.d0+2.d0*dble(j)/dble(10)
!!      xi(3)=-1.d0+2.d0*dble(k)/dble(10)
!!      x=fbem_fem_degbeam3d_x(fbem_line4,x_md,v_md,t_md,xi)
!!      write(63,'(6e25.16)') xi, x
!!    end do
!!    end do
!!    end do
!!    !
!!    call fbem_fem_degbeam3d_K(fbem_line4,x_md,v_md,t_md,Dp,4,4,3,Ks)
!!    call fbem_fem_degbeam3d_Q_body(fbem_line4,x_md,v_md,t_md,4,4,Qs)
!!    do i=1,24
!!      write(99,'(24e25.16)') Ks(:,i)
!!      write(98,'(24e25.16)') Qs(:,i)
!!    end do
!    ! Ejecutar en octave:
!!K=load("fort.99");
!!K=K(7:24,7:24);
!!Q=load("fort.98");
!!Q=0.036*Q(7:24,7:24);
!!b=zeros(18,1);
!!b(3)=1;
!!x=linsolve(K,b);
!!x(3)
!!sqrt(eig(K,Q))
!!w=0.01:0.01:6;
!!for i=1:length(w)
!!x=linsolve(K-w(i)**2*Q,b);
!!xt(i)=x(3);
!!end
!!plot(w,abs(xt))
!!ylim([0 10.0])
!!    deflexion: 0.115 inch
!!    primera frecuencia natural: 5.6 ¿unidades, rad/s?
!  end subroutine test_3dbeamdeg

end module fbem_fem_beams
