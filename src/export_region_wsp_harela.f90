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

subroutine export_region_wsp_harela(kf,kr,c1,c2)

  ! Fortran 2003 intrinsic module
  use iso_fortran_env

  ! fbem modules
  use fbem_data_structures
  use fbem_string_handling
  use fbem_numerical

  use csv_module

  ! Module of problem variables
  use problem_variables

  ! No implicit variables
  implicit none

  ! I/O variables
  integer                           :: kf
  integer                           :: kr
  complex(kind=real64)              :: c1
  complex(kind=real64)              :: c2
  ! Local
  real(kind=real64)                       :: omega
  character(len=fbem_filename_max_length) :: tmp_filename
  type(csv_file)                          :: file_csv
  logical                                 :: status_ok, do_append
  character(len=3)                        :: complex_str1, complex_str2
  character(len=5)                        :: freq_str

  ! Frequency
  omega=frequency(kf)

  call file_csv%initialize
  write(tmp_filename,*) trim(output_filename),'.wsp.','region.',region(kr)%id,'.csv'
  call fbem_trimall(tmp_filename)
  do_append=.true.
  if (kf.eq.1) do_append=.false.
  call file_csv%open(trim(tmp_filename),n_cols=6,status_ok=status_ok,append=do_append)
  if (status_ok) then
    if (kf.eq.1) then
      if (frequency_units.eq.'f') then
        freq_str='Hz'
      else
        freq_str='rad/s'
      end if
      if (complex_notation.eq.2) then
        complex_str1='Re'
        complex_str2='Im'
      else
        complex_str1='Abs'
        complex_str2='Arg'
      end if
      call file_csv%add(['Frequency index','Frequency value ['//trim(freq_str)//']',&
                         trim(complex_str1)//'(c1)',trim(complex_str2)//'(c1)',&
                         trim(complex_str1)//'(c2)',trim(complex_str2)//'(c2)'],trim_str=.true.)
      call file_csv%next_row()
    end if
    call file_csv%add(kf)
    if (frequency_units.eq.'f') then
      call file_csv%add(omega*c_1_2pi)
    else
      call file_csv%add(omega)
    end if
    if (complex_notation.eq.2) then
      call file_csv%add(dreal(c1),real_fmt='('//fmt_real//')')
      call file_csv%add(dimag(c1),real_fmt='('//fmt_real//')')
      call file_csv%add(dreal(c2),real_fmt='('//fmt_real//')')
      call file_csv%add(dimag(c2),real_fmt='('//fmt_real//')')
    else
      call file_csv%add(      abs(c1),real_fmt='('//fmt_real//')')
      call file_csv%add(fbem_zarg(c1),real_fmt='('//fmt_real//')')
      call file_csv%add(      abs(c2),real_fmt='('//fmt_real//')')
      call file_csv%add(fbem_zarg(c2),real_fmt='('//fmt_real//')')
    end if
    call file_csv%next_row()
  else
    call fbem_error_message(error_unit,0,trim(tmp_filename),0,'error when opening the file.')
  end if
  call file_csv%close(status_ok)
  if (.not.status_ok) call fbem_error_message(error_unit,0,trim(tmp_filename),0,'error when closing the file.')

end subroutine export_region_wsp_harela
