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

subroutine help

  ! Fortran 2003 intrinsic module
  use iso_fortran_env

  ! fbem modules
  use fbem_string_handling

  ! String for write format
  character(len=80) :: fmtstr
  character(len=80) :: ver

  fmtstr='(a80)'
  ! Terminal columns
  !                                   1         2         3         4         5         6         7         8
  !                          12345678901234567890123456789012345678901234567890123456789012345678901234567890
  write(output_unit,fmtstr) 'Usage:                                                                          '
  write(output_unit,fmtstr) '  multifebe [options]                                                           '
  write(output_unit,*)
  write(output_unit,fmtstr) 'Options:                                                                        '
  write(output_unit,fmtstr) '  -i, --input STRING        Input file name           [required]                '
  write(output_unit,fmtstr) '  -o, --output STRING       Output files basename     [default: input file name]'
  write(output_unit,fmtstr) '  -m, --memory INTEGER      Max. memory allowed (GB)  [default: 0 (unlimited)]  '
  write(output_unit,fmtstr) '  -b, --verbose INTEGER     Verbose level: 0-10       [default: 1]              '
  write(output_unit,fmtstr) '  -v, --version             Program version                                     '
  write(output_unit,fmtstr) '  -h, --help                This help                                           '
  write(output_unit,*)

  write(output_unit,'(a17)') 'Compilation info:'
#ifdef __GFORTRAN__
  ver=__VERSION__
  call fbem_trim(ver)
  write(fmtstr,*) '(a12,a',len_trim(ver),',a3,1x,a11,1x,a8)'
  call fbem_trimall(fmtstr)
  write(output_unit,fmtstr) 'GNU Fortran ',trim(ver),' on',__DATE__, __TIME__
#endif
#ifdef __INTEL_COMPILER
  ver=__INTEL_COMPILER
  call fbem_trim(ver)
  write(fmtstr,*) '(a14,a',len_trim(ver),',a3,1x,a11,1x,a8)'
  call fbem_trimall(fmtstr)
  write(output_unit,fmtstr) 'Intel Fortran ',trim(ver),' on',__DATE__, __TIME__
#endif
  write(output_unit,*)
  stop

end subroutine help
