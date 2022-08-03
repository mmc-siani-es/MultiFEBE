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
!! <b> This module implements routines for exporting gmsh files. </b>
module fbem_gmsh

  ! Fortran 2003 intrinsic module
  use iso_fortran_env
  ! fbem modules
  use fbem_string_handling

  ! No implicit variables are allowed in the module
  implicit none

  ! By default all is private
  private

  public :: fbem_export_gmsh_fmt_real
  public :: fbem_export_gmsh_NodeData
  public :: fbem_export_gmsh_ElementData
  public :: fbem_export_gmsh_ElementNodeData

  ! Exporting precision
  character(len=fbem_fmtstr) :: fmt_real='e16.8e2'

contains

  subroutine fbem_export_gmsh_fmt_real(fmt_real_val)
    implicit none
    character(len=*) :: fmt_real_val
    fmt_real=fmt_real_val
  end subroutine

  subroutine fbem_export_gmsh_NodeData(fileunit,viewname,timeval,timestep,ncomp,nmax,nnode,nodeid,nodeval)
    implicit none
    ! I/O
    integer           :: fileunit
    character(len=*)  :: viewname
    real(kind=real64) :: timeval
    integer           :: timestep
    integer           :: ncomp
    integer           :: nmax
    integer           :: nnode
    integer           :: nodeid(nmax)
    real(kind=real64) :: nodeval(9,nmax)
    ! Local
    integer                    :: k
    character(len=fbem_fmtstr) :: fmt1, fmt2, fmt_integer
    !
    ! From gmsh documentation:
    !
    !  $NodeData
    !    numStringTags(ASCII int)
    !    stringTag(string) ...
    !    numRealTags(ASCII int)
    !    realTag(ASCII double) ...
    !    numIntegerTags(ASCII int)
    !    integerTag(ASCII int) ...
    !    nodeTag(size_t) value(double) ...
    !    ...
    !  $EndNodeData
    !
    ! check
    if ((ncomp.ne.1).and.(ncomp.ne.3).and.(ncomp.ne.9)) then
      call fbem_error_message(error_unit,0,__FILE__,__LINE__,'invalid value of ncomp')
    end if
    write(fmt_integer,*) 'i', fbem_nchar_int(maxval(nodeid))+1
    call fbem_trimall(fmt_integer)
    ! write
    write(fileunit,'(a9)') '$NodeData'
    write(fileunit,'(a1)' ) '1'
    write(fmt1,*) '(a1,a',len_trim(viewname),'a1)'
    call fbem_trimall(fmt1)
    write(fileunit,fmt1) '"',viewname,'"'
    write(fileunit,'(a1)' ) '1'
    write(fmt1,*) '(',trim(fmt_real),')'
    call fbem_trimall(fmt1)
    write(fileunit,fmt1) timeval
    write(fileunit,'(a1)') '3'
    write(fmt1,*) '(i',fbem_nchar_int(timestep),')'
    call fbem_trimall(fmt1)
    write(fileunit,fmt1) timestep
    write(fileunit,'(i1)') ncomp
    write(fmt1,*) '(i',fbem_nchar_int(nnode),')'
    call fbem_trimall(fmt1)
    write(fileunit,fmt1) nnode
    do k=1,nnode
      write(fmt1,*) '(i',fbem_nchar_int(nodeid(k)),',',ncomp,trim(fmt_real),')'
      call fbem_trimall(fmt1)
      write(fileunit,fmt1) nodeid(k), nodeval(1:ncomp,k)
    end do
    write(fileunit,'(a12)') '$EndNodeData'
  end subroutine fbem_export_gmsh_NodeData

  subroutine fbem_export_gmsh_ElementData(fileunit,viewname,timeval,timestep,ncomp,nmax,nelement,elementid,elementval)
    implicit none
    ! I/O
    integer           :: fileunit
    character(len=*)  :: viewname
    real(kind=real64) :: timeval
    integer           :: timestep
    integer           :: ncomp
    integer           :: nmax
    integer           :: nelement
    integer           :: elementid(nmax)
    real(kind=real64) :: elementval(9,nmax)
    ! Local
    integer                    :: k
    character(len=fbem_fmtstr) :: fmt1, fmt2, fmt_integer
    !
    ! From gmsh documentation:
    !
    !  $ElementData
    !    numStringTags(ASCII int)
    !    stringTag(string) ...
    !    numRealTags(ASCII int)
    !    realTag(ASCII double) ...
    !    numIntegerTags(ASCII int)
    !    integerTag(ASCII int) ...
    !    elementTag(size_t) value(double) ...
    !    ...
    !  $EndElementData
    !
    ! check
    if ((ncomp.ne.1).and.(ncomp.ne.3).and.(ncomp.ne.9)) then
      call fbem_error_message(error_unit,0,__FILE__,__LINE__,'invalid value of ncomp')
    end if
    write(fmt_integer,*) 'i', fbem_nchar_int(maxval(elementid))+1
    call fbem_trim2b(fmt_integer)
    ! write
    write(fileunit,'(a12)') '$ElementData'
    write(fileunit,'(a1)' ) '1'
    write(fmt1,*) '(a1,a',len_trim(viewname),'a1)'
    call fbem_trimall(fmt1)
    write(fileunit,fmt1) '"',viewname,'"'
    write(fileunit,'(a1)' ) '1'
    write(fmt1,*) '(',trim(fmt_real),')'
    call fbem_trimall(fmt1)
    write(fileunit,fmt1) timeval
    write(fileunit,'(a1)') '3'
    write(fmt1,*) '(i',fbem_nchar_int(timestep),')'
    call fbem_trimall(fmt1)
    write(fileunit,fmt1) timestep
    write(fileunit,'(i1)') ncomp
    write(fmt1,*) '(i',fbem_nchar_int(nelement),')'
    call fbem_trimall(fmt1)
    write(fileunit,fmt1) nelement
    do k=1,nelement
      write(fmt1,*) '(i',fbem_nchar_int(elementid(k)),',',ncomp,trim(fmt_real),')'
      call fbem_trimall(fmt1)
      write(fileunit,fmt1) elementid(k), elementval(1:ncomp,k)
    end do
    write(fileunit,'(a15)') '$EndElementData'
  end subroutine fbem_export_gmsh_ElementData

  subroutine fbem_export_gmsh_ElementNodeData(fileunit,viewname,timeval,timestep,ncomp,nemax,nnemax,nelem,elemid,nnodeselem,elemval)
    implicit none
    ! I/O
    integer           :: fileunit
    character(len=*)  :: viewname
    real(kind=real64) :: timeval
    integer           :: timestep
    integer           :: ncomp
    integer           :: nemax
    integer           :: nnemax
    integer           :: nelem
    integer           :: elemid(nemax)
    integer           :: nnodeselem(nemax)
    real(kind=real64) :: elemval(9,nnemax,nemax)
    ! Local
    integer                    :: k
    character(len=fbem_fmtstr) :: fmt1, fmt2, fmt_integer
    !
    ! From gmsh documentation:
    !
    ! $ElementNodeData
    !   numStringTags(ASCII int)
    !   stringTag(string) ...
    !   numRealTags(ASCII int)
    !   realTag(ASCII double) ...
    !   numIntegerTags(ASCII int)
    !   integerTag(ASCII int) ...
    !   elementTag(size_t) numNodesPerElement(int) value(double) ...
    !   ...
    ! $EndElementNodeData
    !
    ! check
    if ((ncomp.ne.1).and.(ncomp.ne.3).and.(ncomp.ne.9)) then
      call fbem_error_message(error_unit,0,__FILE__,__LINE__,'invalid value of ncomp')
    end if
    write(fmt_integer,*) 'i', fbem_nchar_int(maxval(elemid))+1
    call fbem_trim2b(fmt_integer)
    ! write
    write(fileunit,'(a16)') '$ElementNodeData'
    write(fileunit,'(a1)' ) '1'
    write(fmt1,*) '(a1,a',len_trim(viewname),'a1)'
    call fbem_trimall(fmt1)
    write(fileunit,fmt1) '"',viewname,'"'
    write(fileunit,'(a1)' ) '1'
    write(fmt1,*) '(',trim(fmt_real),')'
    call fbem_trimall(fmt1)
    write(fileunit,fmt1) timeval
    write(fileunit,'(a1)') '3'
    write(fmt1,*) '(i',fbem_nchar_int(timestep),')'
    call fbem_trimall(fmt1)
    write(fileunit,fmt1) timestep
    write(fileunit,'(i2)') ncomp
    write(fmt1,*) '(i',fbem_nchar_int(nelem),')'
    call fbem_trimall(fmt1)
    write(fileunit,fmt1) nelem
    do k=1,nelem
      write(fmt1,*) '(i',fbem_nchar_int(elemid(k)),',1x,i',fbem_nchar_int(nnodeselem(k)),',',ncomp*nnodeselem(k),trim(fmt_real),')'
      call fbem_trimall(fmt1)
      write(fileunit,fmt1) elemid(k), nnodeselem(k), elemval(1:ncomp,1:nnodeselem(k),k)
    end do
    write(fileunit,'(a19)') '$EndElementNodeData'
  end subroutine fbem_export_gmsh_ElementNodeData

end module fbem_gmsh
