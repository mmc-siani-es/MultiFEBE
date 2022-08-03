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

!! This subroutine solves a Linear System of Equations (LSE) using several subroutines from the Lapack library. It is possible to
!! perform the scaling of the LSE to improve the condition number, estimate the condition number and/or refine the solution. For
!! complex numbers of <tt>kind=real64</tt> and multiple right hand sides.
subroutine solve_lse_c(n_dof,A,ipiv,Aodim,Ao,equed,r,c,n_rhs,b,factorize,scaling,condition,refine)
  ! Fortran 2003 intrinsic module
  use iso_fortran_env
  ! fbem modules
  use fbem_string_handling
  ! No implicit variables are allowed
  implicit none
  ! I/O
  integer                           :: n_dof           !! Number of degrees of freedom
  complex(kind=real64)              :: A(n_dof,n_dof)  !! Matrix A
  integer                           :: ipiv(n_dof)     !! Row order
  integer                           :: Aodim           !! If condition||refine==T, Aodim should be n_dof, otherwise Ao not used.
  complex(kind=real64)              :: Ao(Aodim,Aodim) !! Matrix A unfactorized
  character(len=1)                  :: equed           !! Scaling applied
  real(kind=real64)                 :: r(n_dof)        !! Row scaling factors
  real(kind=real64)                 :: c(n_dof)        !! Column scaling factors
  integer                           :: n_rhs           !! Number of right hand sides
  complex(kind=real64)              :: b(n_dof,n_rhs)  !! Vector b
  logical                           :: factorize       !! True if A has to be factorized
  logical                           :: scaling         !! Perform scaling
  logical                           :: condition       !! Estimate the condition number
  logical                           :: refine          !! Perform solution refinement
  ! Local
  complex(kind=real64), allocatable :: bcopy(:,:)      !! Copy of vector b required by the refine stage
  real(kind=real64)                 :: rowcnd
  real(kind=real64)                 :: colcnd
  real(kind=real64)                 :: amax
  real(kind=real64)                 :: summ
  character(len=1)                  :: norm
  real(kind=real64)                 :: anorm
  real(kind=real64)                 :: rcond
  real(kind=real64)                 :: rwork(n_dof)
  complex(kind=real64)              :: work(4*n_dof)
  integer                           :: info
  character(len=1)                  :: trans
  integer                           :: i, j
  real(kind=real64)                 :: ferr(n_rhs)
  real(kind=real64)                 :: berr(n_rhs)

  if (verbose_level.ge.1)  call fbem_timestamp_w_message(output_unit,2,'START solving LSE')

  ! Initialize Ao if required
  if (condition.or.refine) then
    if (Aodim.eq.n_dof) then
      if (factorize) then
        Ao=A
      end if
    else
      call fbem_error_message(error_unit,0,__FILE__,__LINE__,'invalid value of Aodim')
    end if
  end if

  ! ------------- !
  ! SCALING STAGE !
  ! ------------- !

  if (scaling) then
    ! Perform scaling of A (and Ao) only if A is not already factorized
    if (factorize) then
      call zgeequ(n_dof,n_dof,A,n_dof,r,c,rowcnd,colcnd,amax,info)
      if (info.eq.0) then
        if (verbose_level.ge.2) then
          call fbem_timestamp_w_message(output_unit,2,'The scaling factors for the A matrix has been obtained correctly')
        end if
      else
        call fbem_timestamp_message(output_unit,2)
        write(output_unit,'(a13,i11)') 'zgeequ info: ', info
        call fbem_error_message(error_unit,0,__FILE__,__LINE__,'scaling has failed')
      end if
      call zlaqge(n_dof,n_dof,A,n_dof,r,c,rowcnd,colcnd,amax,equed)
      if (condition.or.refine) Ao=A
      select case (equed)
        case ('N')
          if (verbose_level.ge.2) call fbem_timestamp_w_message(output_unit,2,'No equilibration')
        case ('R')
          if (verbose_level.ge.2) call fbem_timestamp_w_message(output_unit,2,'Row equilibration')
        case ('C')
          if (verbose_level.ge.2) call fbem_timestamp_w_message(output_unit,2,'Column equilibration')
        case ('B')
          if (verbose_level.ge.2) call fbem_timestamp_w_message(output_unit,2,'Both row and column equilibration')
      end select
    else
      ! Check equed from input
      if ((equed.ne.'N').and.(equed.ne.'R').and.(equed.ne.'C').and.(equed.ne.'B')) then
        call fbem_error_message(error_unit,0,__FILE__,__LINE__,'invalid value of equed')
      end if
    end if
    if ((equed.eq.'R').or.(equed.eq.'B')) then
      do i=1,n_rhs
        b(:,i)=r(:)*b(:,i)
      end do
    end if
  end if

  ! ------------------- !
  ! FACTORIZATION STAGE !
  ! ------------------- !

  if (factorize) then
    call zgetrf(n_dof,n_dof,A,n_dof,ipiv,info)
    if (info.eq.0) then
      if (verbose_level.ge.2) then
        call fbem_timestamp_w_message(output_unit,2,'Linear system of equations LU-factorized')
      end if
    else
      call fbem_timestamp_message(output_unit,2)
      write(output_unit,'(a13,i11)') 'zgetrf info: ', info
      call fbem_error_message(error_unit,0,__FILE__,__LINE__,'factorization has failed')
    end if
  end if

  ! ---------------------------------- !
  ! CONDITION NUMBER CALCULATION STAGE !
  ! ---------------------------------- !

  if (condition.and.factorize) then
    norm='1'
    !anorm=zlange(norm,n_dof,n_dof,Ao,n_dof,work)
    anorm=0.0d0
    do j=1,n_dof
      summ=0.0d0
      do i=1,n_dof
        summ=summ+abs(Ao(i,j))
      end do
      anorm=max(anorm,summ)
    end do
    call zgecon(norm,n_dof,A,n_dof,anorm,rcond,work,rwork,info)
    if (info.eq.0) then
      if (verbose_level.ge.2) then
        call fbem_timestamp_w_message(output_unit,2,'Condition number correctly estimated')
      end if
    else
      call fbem_timestamp_message(output_unit,2)
      write(output_unit,'(a13,i11)') 'zgecon info: ', info
      call fbem_error_message(error_unit,0,__FILE__,__LINE__,'condition number calculation has failed')
    end if
    if (verbose_level.ge.2) then
      call fbem_timestamp_message(output_unit,2)
      write(output_unit,'(a17,e25.16)') 'Condition number:', 1.0d0/rcond
    end if
  end if

  ! ----------- !
  ! SOLVE STAGE !
  ! ----------- !

  if (refine) then
    allocate (bcopy(n_dof,n_rhs))
    bcopy=b
  end if
  trans='N'
  call zgetrs(trans,n_dof,n_rhs,A,n_dof,ipiv,b,n_dof,info)
  if (info.eq.0) then
    if (verbose_level.ge.2) then
      call fbem_timestamp_w_message(output_unit,2,'Linear system of equations solved')
    end if
  else
    call fbem_timestamp_message(output_unit,2)
    write(output_unit,'(a13,i11)') 'zgetrs info: ', info
    call fbem_error_message(error_unit,0,__FILE__,__LINE__,'solver has failed')
  end if

  ! ------------ !
  ! REFINE STAGE !
  ! ------------ !

  if (refine) then
    call zgerfs(trans,n_dof,n_rhs,Ao,n_dof,A,n_dof,ipiv,bcopy,n_dof,b,n_dof,ferr,berr,work,rwork,info)
    if (info.eq.0) then
      if (verbose_level.ge.2) then
        do i=1,n_rhs
          call fbem_timestamp_message(output_unit,2)
          write(output_unit,'(a30,i11,a5,e25.16,a5,e25.16)') 'Solution refinement done: rhs=',i,', berr=',berr(i),', ferr=',ferr(i)
        end do
      end if
    else
      call fbem_timestamp_message(output_unit,2)
      write(output_unit,'(a13,i11)') 'zgerfs info: ', info
      call fbem_error_message(error_unit,0,__FILE__,__LINE__,'solution refinement has failed')
    end if
    deallocate (bcopy)
  end if

  ! Unscaling of b
  if (scaling) then
    if ((equed.eq.'C').or.(equed.eq.'B')) then
      do i=1,n_rhs
        b(:,i)=c(:)*b(:,i)
      end do
    end if
  end if

  if (verbose_level.ge.1) call fbem_timestamp_w_message(output_unit,2,'END solving LSE')

end subroutine solve_lse_c
