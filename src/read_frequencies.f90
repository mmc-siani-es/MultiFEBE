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
!! <b> Subroutine that reads the frequencies from a file and build
!!     the vector of frequencies. </b>
!!
!!>
!! 'Hz' or 'rad/s'
!! 'list', 'lin', 'log10', 'sound_octave', 'sound_1/3-octave'
!! n_frequencies
!! frequency(1) (minimum frequency for mode 'lin' and 'log10')
!! ... (if mode is 'list')
!! frequency(n_frequencies) (maximum frequency for mode 'lin' and 'log10')
!!<
subroutine read_frequencies(fileunit)

  ! fbem module
  use fbem_numerical
  use fbem_string_handling
  use fbem_acoustics

  ! Problem variables module
  use problem_variables

  implicit none
  ! I/O
  integer                                 :: fileunit
  ! Local
  character(len=fbem_stdcharlen)          :: section_name      ! Name of the section
  character(len=fbem_file_record_length)  :: line
  integer                                 :: i, nf
  integer                                 :: specification_mode
  real(kind=real64)                       :: delta, minf, maxf
  logical                                 :: found

  section_name='frequencies'
  if (verbose_level.ge.2) call fbem_timestamp_w_message(output_unit,2,'SEARCHING section ['//trim(section_name)//']')
  call fbem_search_section(fileunit,section_name,found)
  if (found) then

    if (verbose_level.ge.2) call fbem_timestamp_w_message(output_unit,2,'START reading section ['//trim(section_name)//']')

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
    if (trim(line).eq.'list'            ) specification_mode=1
    if (trim(line).eq.'lin'             ) specification_mode=2
    if (trim(line).eq.'log10'           ) specification_mode=3
    if (trim(line).eq.'log'             ) specification_mode=3
    if (trim(line).eq.'sound_octave'    ) specification_mode=4
    if (trim(line).eq.'sound_1/3-octave') specification_mode=5
    if (specification_mode.eq.0) then
      call fbem_error_message(error_unit,0,'specification_mode',0,'the frequency specification mode is not valid')
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
      !
      ! sound_octave
      !
      case (4)
        ! Read the number of frequencies
        read(fileunit,*) n_frequencies
        ! Read minimum and maximum frequencies
        read(fileunit,*) minf
        read(fileunit,*) maxf
        ! If n_frequencies == 0 (all octave frequencies between minf and maxf are analyzed)
        if (n_frequencies.eq.0) then
          n_frequencies = 0
          do i=1,size(fbem_octave_f)
            if ((fbem_octave_f(i).ge.minf).and.(fbem_octave_f(i).le.maxf)) then
              n_frequencies = n_frequencies + 1
            end if
          end do
          allocate(frequency(n_frequencies))
          n_frequencies = 0
          do i=1,size(fbem_octave_f)
            if ((fbem_octave_f(i).ge.minf).and.(fbem_octave_f(i).le.maxf)) then
              n_frequencies = n_frequencies + 1
              frequency(n_frequencies) = fbem_octave_f(i)
            end if
          end do
        ! If n_frequencies > 0 (up to n_frequencies between minf and maxf are analyzed starting from minf)
        elseif (n_frequencies.gt.0) then
          nf = 0
          do i=1,size(fbem_octave_f)
            if ((fbem_octave_f(i).ge.minf).and.(fbem_octave_f(i).le.maxf)) then
              nf = nf + 1
            end if
          end do
          n_frequencies = min(n_frequencies,nf)
          allocate(frequency(n_frequencies))
          n_frequencies = 0
          do i=1,size(fbem_octave_f)
            if ((fbem_octave_f(i).ge.minf).and.(fbem_octave_f(i).le.maxf)) then
              n_frequencies = n_frequencies + 1
              frequency(n_frequencies) = fbem_octave_f(i)
            end if
          end do
        ! If n_frequencies > 0 (up to n_frequencies between minf and maxf are analyzed starting from maxf)
        elseif (n_frequencies.lt.0) then
          nf = 0
          do i=size(fbem_octave_f),1
            if ((fbem_octave_f(i).ge.minf).and.(fbem_octave_f(i).le.maxf)) then
              nf = nf + 1
            end if
          end do
          n_frequencies = min(n_frequencies,nf)
          allocate(frequency(n_frequencies))
          n_frequencies = 0
          do i=size(fbem_octave_f),1
            if ((fbem_octave_f(i).ge.minf).and.(fbem_octave_f(i).le.maxf)) then
              n_frequencies = n_frequencies + 1
              frequency(n_frequencies) = fbem_octave_f(i)
            end if
          end do

        end if
      !
      ! sound_1/3-octave
      !
      case (5)
        ! Read the number of frequencies
        read(fileunit,*) n_frequencies
        ! Read minimum and maximum frequencies
        read(fileunit,*) minf
        read(fileunit,*) maxf
        ! If n_frequencies == 0 (all octave frequencies between minf and maxf are analyzed)
        if (n_frequencies.eq.0) then
          n_frequencies = 0
          do i=1,size(fbem_onethird_octave_f)
            if ((fbem_onethird_octave_f(i).ge.minf).and.(fbem_onethird_octave_f(i).le.maxf)) then
              n_frequencies = n_frequencies + 1
            end if
          end do
          allocate(frequency(n_frequencies))
          n_frequencies = 0
          do i=1,size(fbem_onethird_octave_f)
            if ((fbem_onethird_octave_f(i).ge.minf).and.(fbem_onethird_octave_f(i).le.maxf)) then
              n_frequencies = n_frequencies + 1
              frequency(n_frequencies) = fbem_onethird_octave_f(i)
            end if
          end do
        ! If n_frequencies > 0 (up to n_frequencies between minf and maxf are analyzed starting from minf)
        elseif (n_frequencies.gt.0) then
          nf = 0
          do i=1,size(fbem_onethird_octave_f)
            if ((fbem_onethird_octave_f(i).ge.minf).and.(fbem_onethird_octave_f(i).le.maxf)) then
              nf = nf + 1
            end if
          end do
          n_frequencies = min(n_frequencies,nf)
          allocate(frequency(n_frequencies))
          n_frequencies = 0
          do i=1,size(fbem_onethird_octave_f)
            if ((fbem_onethird_octave_f(i).ge.minf).and.(fbem_onethird_octave_f(i).le.maxf)) then
              n_frequencies = n_frequencies + 1
              frequency(n_frequencies) = fbem_onethird_octave_f(i)
            end if
          end do
        ! If n_frequencies > 0 (up to n_frequencies between minf and maxf are analyzed starting from maxf)
        elseif (n_frequencies.lt.0) then
          nf = 0
          do i=size(fbem_onethird_octave_f),1
            if ((fbem_onethird_octave_f(i).ge.minf).and.(fbem_onethird_octave_f(i).le.maxf)) then
              nf = nf + 1
            end if
          end do
          n_frequencies = min(n_frequencies,nf)
          allocate(frequency(n_frequencies))
          n_frequencies = 0
          do i=size(fbem_onethird_octave_f),1
            if ((fbem_onethird_octave_f(i).ge.minf).and.(fbem_onethird_octave_f(i).le.maxf)) then
              n_frequencies = n_frequencies + 1
              frequency(n_frequencies) = fbem_onethird_octave_f(i)
            end if
          end do
        end if
    end select

    ! If input/output units are Hz, then the vector of frequencies must be converted to rad/s.
    if (frequency_units.eq.'f') then
      do i=1,n_frequencies
        frequency(i)=c_2pi*frequency(i)
      end do
    end if

    if (verbose_level.ge.2) call fbem_timestamp_w_message(output_unit,2,'END reading section ['//trim(section_name)//']')

  else

    call fbem_error_message(error_unit,0,trim(section_name),0,'this section is required')

  end if

end subroutine read_frequencies
