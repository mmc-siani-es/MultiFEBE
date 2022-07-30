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


!! @author Jacob David Rodriguez Bordon (jacobdavid.rodriguezbordon@ulpgc.es)
!!
!! @version 2.0
!!
!! <b> Subroutine that reads the frequencies from a file. </b>
!!
!!>
!! 'Hz' or 'rad/s'
!! 'list' or 'lin' or 'log10'
!! n_frequencies
!! frequency(1) (minimum frequency for mode 'lin' and 'log10')
!! ... (if mode is 'list')
!! frequency(n_frequencies) (maximum frequency for mode 'lin' and 'log10')
!!<
subroutine read_frequencies(fileunit)

  ! fbem module
  use fbem_numerical
  use fbem_string_handling

  ! Problem variables module
  use problem_variables

  implicit none
  ! I/O
  integer                                 :: fileunit
  ! Local
  character(len=fbem_file_record_length)  :: line
  integer                                 :: i
  integer                                 :: specification_mode
  real(kind=real64)                       :: delta
  logical                                 :: found

  if (verbose_level.ge.2) call fbem_timestamp_w_message(output_unit,2,'SEARCHING section [frequencies]')
  call fbem_search_section(fileunit,'frequencies',found)
  if (found) then

    if (verbose_level.ge.2) call fbem_timestamp_w_message(output_unit,2,'START reading section [frequencies]')

    ! Read units
    read(fileunit,'(a)') line
    call fbem_trim2b(line)
    frequency_units=''
    if (trim(line).eq.'Hz'   ) frequency_units='f'
    if (trim(line).eq.'rad/s') frequency_units='w'
    if (trim(frequency_units).eq.'') then
      call fbem_error_message(error_unit,0,'frequency_units',0,'the frequency units can be only "Hz" or "rad/s"')
    end if

    ! Read the frequency specification mode
    read(fileunit,'(a)') line
    call fbem_trim2b(line)
    specification_mode=0
    if (trim(line).eq.'list' ) specification_mode=1
    if (trim(line).eq.'lin'  ) specification_mode=2
    if (trim(line).eq.'log10') specification_mode=3
    if (trim(line).eq.'log'  ) specification_mode=3
    if (specification_mode.eq.0) then
      call fbem_error_message(error_unit,0,'specification_mode',0,'the frequency specification mode is "list", "lin" or "log"')
    end if

    ! Switch between the frequency specification modes
    select case (specification_mode)
      !
      ! list
      !
      case (1)
        ! Read the number of frequencies
        read(fileunit,*) n_frequencies
        if (n_frequencies.le.0) then
          call fbem_error_message(error_unit,0,'n_frequencies',0,'the number of frequencies must be >=1')
        end if
        ! Read frequencies
        allocate(frequency(n_frequencies))
        do i=1,n_frequencies
          read(fileunit,*) frequency(i)
        end do
      !
      ! lin
      !
      case (2)
        ! Read the number of frequencies
        read(fileunit,*) n_frequencies
        if (n_frequencies.lt.2) then
          call fbem_error_message(error_unit,0,'n_frequencies',0,'the number of frequencies must be >=2')
        end if
        ! Read minimum and maximum frequencies
        allocate(frequency(n_frequencies))
        read(fileunit,*) frequency(1)
        read(fileunit,*) frequency(n_frequencies)
        if (frequency(1).ge.frequency(n_frequencies)) then
          call fbem_error_message(error_unit,0,'frequency',0,'the range of frequencies is invalid')
        end if
        ! Build the vector of frequencies
        delta=(frequency(n_frequencies)-frequency(1))/dble(n_frequencies-1)
        do i=2,n_frequencies-1
          frequency(i)=frequency(1)+delta*dble(i-1)
        end do
      !
      ! log
      !
      case (3)
        ! Read the number of frequencies
        read(fileunit,*) n_frequencies
        if (n_frequencies.lt.2) then
          call fbem_error_message(error_unit,0,'n_frequencies',0,'the number of frequencies must be >=2')
        end if
        ! Read minimum and maximum frequencies
        allocate(frequency(n_frequencies))
        read(fileunit,*) frequency(1)
        read(fileunit,*) frequency(n_frequencies)
        if (frequency(1).ge.frequency(n_frequencies)) then
          call fbem_error_message(error_unit,0,'frequency',0,'the range of frequencies is invalid')
        end if
        ! Build the vector of frequencies
        delta=(dlog10(frequency(n_frequencies))-dlog10(frequency(1)))/dble(n_frequencies-1)
        do i=2,n_frequencies-1
          frequency(i)=10.0d0**(dlog10(frequency(1))+delta*dble(i-1))
        end do
    end select

    ! If input units are Hz, then the vector of frequencies must be converted to rad/s.
    if (frequency_units.eq.'f') then
      do i=1,n_frequencies
        frequency(i)=c_2pi*frequency(i)
      end do
    end if

    if (verbose_level.ge.2) call fbem_timestamp_w_message(output_unit,2,'END reading section [frequencies]')

  else

    call fbem_error_message(error_unit,0,'[frequencies]',0,'this section is required')

  end if

end subroutine read_frequencies
