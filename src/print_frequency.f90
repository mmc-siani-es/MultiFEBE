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

subroutine print_frequency(selected_unit,print_mode,kf)

  use iso_fortran_env
  use fbem_string_handling
  use problem_variables

  implicit none

  integer                    :: selected_unit
  integer                    :: print_mode
  integer                    :: kf

  character(len=fbem_fmtstr) :: fmtstr

  select case (print_mode)

    case (1)

      if (verbose_level.ge.1)  then
        call fbem_timestamp_message(selected_unit,2)
        write(fmtstr,*) '(a,i',fbem_nchar_int(n_frequencies),',a,i',fbem_nchar_int(n_frequencies),',a,f5.1,a)'
        call fbem_trimall(fmtstr)
        write(selected_unit,fmtstr) 'Frequency ',kf,'/',n_frequencies, ' (', dble(kf)/dble(n_frequencies)*100., '%)'
      end if

    case default

      call fbem_error_message(error_unit,0,__FILE__,__LINE__,'print_mode not valid')

  end select

end subroutine print_frequency
