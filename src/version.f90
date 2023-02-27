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

subroutine version

  ! Fortran 2003 intrinsic module
  use iso_fortran_env

  ! String for write format
  character(len=5) :: fmtstr

  fmtstr='(a80)'
  ! Terminal columns
  !                                   1         2         3         4         5         6         7         8
  !                          12345678901234567890123456789012345678901234567890123456789012345678901234567890
  write(output_unit,fmtstr) 'MultiFEBE 2.0.3                                                                 '
  write(output_unit,fmtstr) 'Copyright (C) 2014-2022 Universidad de Las Palmas de Gran Canaria:              '
  write(output_unit,fmtstr) '                        Jacob D.R. Bordon                                       '
  write(output_unit,fmtstr) '                        Guillermo M. Alamo                                      '
  write(output_unit,fmtstr) '                        Juan J. Aznarez                                         '
  write(output_unit,fmtstr) '                        Orlando Maeso                                           '
  write(output_unit,*)
  write(output_unit,fmtstr) 'Description: Multi-domain integrated Finite Element and Boundary Element solver '
  write(output_unit,*)
  write(output_unit,fmtstr) 'Link: https://github.com/mmc-siani-es/MultiFEBE                                 '
  write(output_unit,*)
  write(output_unit,fmtstr) 'License GPLv2: GNU GPL version 2 or later <https://gnu.org/licenses/gpl.html>   '
  write(output_unit,fmtstr) 'This is free software; see the source for copying conditions.  There is NO      '
  write(output_unit,fmtstr) 'warranty; not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.     '
  write(output_unit,*)

end subroutine version
