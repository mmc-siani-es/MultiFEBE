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
!! <b>This module has string handling parameters, subroutines and functions.</b>
module fbem_string_handling

  ! Fortran 2003 intrinsic module
  use iso_fortran_env

  ! No implicit variables are allowed in the module
  implicit none

  ! By default all is private
  private

  ! Public functions

  !! Verbose level (default 1) for warning and error messages subroutines. If <tt>verbose_level=0</tt> nothing is printed.
  integer, public :: verbose_level=1

  ! String routines
  public :: fbem_timestamp_message
  public :: fbem_timestamp_w_message
  public :: fbem_warning_message
  public :: fbem_error_message
  public :: fbem_text_style
  public :: fbem_trim
  public :: fbem_trimall
  public :: fbem_trim2b
  public :: fbem_nchar_int
  public :: fbem_count_words
  public :: fbem_extract_word
  public :: fbem_search_section
  public :: fbem_search_section_gmsh
  public :: fbem_search_keyword
  ! TO-DO: in another module called fbem_file_handling or similar
  public :: fbem_get_valid_unit
  public :: fbem_get_current_working_directory
  public :: fbem_get_os_path_type
  public :: fbem_path_is_absolute
  public :: fbem_path_is_relative
  public :: fbem_get_dirname
  public :: fbem_file_exists

  ! Parameters
  integer, parameter, public :: fbem_stdcharlen         =32   !! Standard length for entities names
  integer, parameter, public :: fbem_string_max_length  =2048 !! Maximum string length for subroutines
  integer, parameter, public :: fbem_filename_max_length=255  !! Maximum file name length
  integer, parameter, public :: fbem_path_max_length    =4096 !! Maximum path length
  integer, parameter, public :: fbem_file_record_length =8192 !! File record length
  integer, parameter, public :: fbem_fmtstr             =1024 !! Format string length
  ! ANSI escape codes
  character(len=1), parameter, public :: fbem_c_esc = achar(27)
  character(len=2), parameter, public :: fbem_c_start = fbem_c_esc // '['
  character(len=1), parameter, public :: fbem_c_end = 'm'
  character(len=*), parameter, public :: fbem_c_clear = fbem_c_start // '0' // fbem_c_end
  character(len=*), parameter, public :: fbem_st_bold = '1'
  character(len=*), parameter, public :: fbem_st_italic = '3'
  character(len=*), parameter, public :: fbem_st_underline = '4'
  character(len=*), parameter, public :: fbem_fg_black = '30'
  character(len=*), parameter, public :: fbem_fg_red = '31'
  character(len=*), parameter, public :: fbem_fg_green = '32'
  character(len=*), parameter, public :: fbem_fg_yellow = '33'
  character(len=*), parameter, public :: fbem_fg_blue = '34'
  character(len=*), parameter, public :: fbem_fg_magenta = '35'
  character(len=*), parameter, public :: fbem_fg_cyan = '36'
  character(len=*), parameter, public :: fbem_fg_white = '37'
  character(len=*), parameter, public :: fbem_bg_black = '40'
  character(len=*), parameter, public :: fbem_bg_red = '41'
  character(len=*), parameter, public :: fbem_bg_green = '42'
  character(len=*), parameter, public :: fbem_bg_yellow = '43'
  character(len=*), parameter, public :: fbem_bg_blue = '44'
  character(len=*), parameter, public :: fbem_bg_magenta = '45'
  character(len=*), parameter, public :: fbem_bg_cyan = '46'
  character(len=*), parameter, public :: fbem_bg_white = '47'

  !! Interface for overloaded functions: <tt>fbem_nchar_int_k8, fbem_nchar_int_k16, fbem_nchar_int_k32, fbem_nchar_int_k64</tt>.
  !! They return the number of characters of a given integer number.
  interface fbem_nchar_int
    module procedure fbem_nchar_int_k8
    module procedure fbem_nchar_int_k16
    module procedure fbem_nchar_int_k32
    module procedure fbem_nchar_int_k64
  end interface

contains

  !! Print colored text. Codes can be separated by semicolor ";"
  function fbem_text_style(str,codes) result(out)
    character(len=*), intent(in) :: str
    character(len=*), intent(in) :: codes
    character(len=:), allocatable :: out
    out = fbem_c_start // codes // fbem_c_end // str // fbem_c_clear
  end function fbem_text_style

  !! Remove beginning and ending blanks of a string
  subroutine fbem_trim(arg)
    implicit none
    character(len=*), intent(inout) :: arg
    arg=trim(adjustl(arg))
  end subroutine fbem_trim

  !! Remove all blanks of a string
  subroutine fbem_trimall(arg)
    implicit none
    character(len=*), intent(inout) :: arg
    integer                         :: i,n,b
    arg=trim(adjustl(arg))
    n=len_trim(arg)
    b=0
    do i=2,n
      if (arg(i:i).eq.' ') then
        b=b+1
      else
        arg((i-b):(i-b))=arg(i:i)
      end if
    end do
    arg=arg(1:n-b)
  end subroutine fbem_trimall

  !! Remove beginning, ending and duplicated blanks of a string
  subroutine fbem_trim2b(arg)
    implicit none
    character(len=*), intent(inout) :: arg
    integer                         :: i,n,b
    logical                         :: pb
    arg=trim(adjustl(arg))
    n=len_trim(arg)
    b=0
    pb=.false.
    do i=2,n
      if (pb.eqv.(.false.)) then
        if (arg(i:i).eq.' ') pb=.true.
        arg((i-b):(i-b))=arg(i:i)
      else
        if (arg(i:i).eq.' ') then
          b=b+1
        else
          pb=.false.
          arg((i-b):(i-b))=arg(i:i)
        end if
      end if
    end do
    arg=arg(1:n-b)
  end subroutine fbem_trim2b

  !! Return the number of characters of a given integer number of kind 8 bits.
  function fbem_nchar_int_k8(vint)
    integer(kind=int8)            :: fbem_nchar_int_k8
    integer(kind=int8), intent(in) :: vint
    if (vint.gt.0) then
      fbem_nchar_int_k8=int(log10(dble(abs(vint))))+1
    else
      if (vint.lt.0) then
        fbem_nchar_int_k8=int(log10(dble(abs(vint))))+2
      else
      if (vint.lt.0) then
        fbem_nchar_int_k8=int(log10(dble(abs(vint))))+2
      else
        fbem_nchar_int_k8=1
      end if
      end if
    end if
  end function fbem_nchar_int_k8

  !! Return the number of characters of a given integer number of kind 16 bits.
  function fbem_nchar_int_k16(vint)
    integer(kind=int16)             :: fbem_nchar_int_k16
    integer(kind=int16), intent(in) :: vint
    if (vint.gt.0) then
      fbem_nchar_int_k16=floor(log10(dble(abs(vint))))+1
    else
      if (vint.lt.0) then
        fbem_nchar_int_k16=floor(log10(dble(abs(vint))))+2
      else
        fbem_nchar_int_k16=1
      end if
    end if
  end function fbem_nchar_int_k16

  !! Return the number of characters of a given integer number of kind 32 bits.
  function fbem_nchar_int_k32(vint)
    integer(kind=int32)             :: fbem_nchar_int_k32
    integer(kind=int32), intent(in) :: vint
    if (vint.gt.0) then
      fbem_nchar_int_k32=floor(log10(dble(vint)))+1
    else
      if (vint.lt.0) then
        fbem_nchar_int_k32=floor(log10(dble(abs(vint))))+2
      else
        fbem_nchar_int_k32=1
      end if
    end if
  end function fbem_nchar_int_k32

  !! Return the number of characters of a given integer number of kind 64 bits.
  function fbem_nchar_int_k64(vint)
    integer(kind=int32)             :: fbem_nchar_int_k64
    integer(kind=int64), intent(in) :: vint
    if (vint.gt.0) then
      fbem_nchar_int_k64=floor(log10(dble(abs(vint))))+1
    else
      if (vint.lt.0) then
        fbem_nchar_int_k64=floor(log10(dble(abs(vint))))+2
      else
        fbem_nchar_int_k64=1
      end if
    end if
  end function fbem_nchar_int_k64

  !! Print a warning message:
  !!
  !! <tt>
  !! Warning (file_name:line_number): message
  !! </tt>
  subroutine fbem_warning_message(selected_unit,verbose_limit,file_name,line_number,message)
    implicit none
    ! I/O
    integer          :: selected_unit !! Unit used to print the message
    integer          :: verbose_limit !! Minimum Verbose level needed to print the message
    character(len=*) :: file_name     !! File name, it is recommended to use the macro <tt>__FILE__</tt>. Also it could be used as a field to print the function/subroutine name.
    integer          :: line_number   !! Line number, it is recommended to use the macro <tt>__LINE__</tt>. Also it could be used as a field to print the warning code.
    character(len=*) :: message       !! Message to print.
    ! Local
    integer                    :: n, m, l, p
    character(len=fbem_fmtstr) :: fmtstr
    character(len=fbem_string_max_length) :: warningstr
    ! Only if verbose_level>=verbose_limit
    if (verbose_level.ge.verbose_limit) then
      write(warningstr,*) fbem_text_style('Warning',fbem_st_bold//';'//fbem_fg_magenta)
      call fbem_trimall(warningstr)
      n=len_trim(file_name)
      m=len_trim(message)
      l=fbem_nchar_int(line_number)
      p=len_trim(warningstr)
      if (line_number.eq.0) l=0
      if ((n+l).gt.0) then
        if ((n.gt.0).and.(l.gt.0)) then
          write(fmtstr,*) '(a',p,',a2,a',n,',a1,i',l,',a3,a',m,')'
          call fbem_trimall(fmtstr)
          write (selected_unit,fmtstr) trim(warningstr),' (',trim(file_name),':',line_number,'): ',trim(message)
        else
          if (n.eq.0) then
            write(fmtstr,*) '(a',p,',a2,i',l,',a3,a',m,')'
            call fbem_trimall(fmtstr)
            write (selected_unit,fmtstr) trim(warningstr),' (',line_number,'): ',trim(message)
          else
            write(fmtstr,*) '(a',p,',a2,a',n,',a3,a',m,')'
            call fbem_trimall(fmtstr)
            write (selected_unit,fmtstr) trim(warningstr),' (',trim(file_name),'): ',trim(message)
          end if
        end if
      else
        write(fmtstr,*) '(a',p,',a2,a',m,')'
        call fbem_trimall(fmtstr)
        write (selected_unit,fmtstr) trim(warningstr),': ',trim(message)
      end if
    end if
  end subroutine fbem_warning_message

  !! Print an error message and stop execution:
  !!
  !! <tt>
  !! Error (file_name:line_number): message
  !! </tt>
  subroutine fbem_error_message(selected_unit,verbose_limit,file_name,line_number,message)
    implicit none
    ! I/O
    integer          :: selected_unit !! Unit used to print the message
    integer          :: verbose_limit !! Minimum verbose level needed to print the message
    character(len=*) :: file_name     !! File name, it is recommended to use the macro <tt>__FILE__</tt>. Also it could be used as a field to print the function/subroutine name.
    integer          :: line_number   !! Line number, it is recommended to use the macro <tt>__LINE__</tt>. Also it could be used as a field to print the warning code.
    character(len=*) :: message       !! Message to print.
    ! Local
    integer                    :: n, m, l, p
    character(len=fbem_fmtstr) :: fmtstr
    character(len=fbem_string_max_length) :: errorstr
    ! Only if verbose_level>=verbose_limit
    if (verbose_level.ge.verbose_limit) then
      write(errorstr,*) fbem_text_style('Error',fbem_st_bold//';'//fbem_fg_red)
      call fbem_trimall(errorstr)
      n=len_trim(file_name)
      m=len_trim(message)
      l=fbem_nchar_int(line_number)
      p=len_trim(errorstr)
      if (line_number.eq.0) l=0
      if ((n+l).gt.0) then
        if ((n.gt.0).and.(l.gt.0)) then
          write(fmtstr,*) '(a',p,',a2,a',n,',a1,i',l,',a3,a',m,')'
          call fbem_trimall(fmtstr)
          write (selected_unit,fmtstr) trim(errorstr),' (',trim(file_name),':',line_number,'): ',trim(message)
        else
          if (n.eq.0) then
            write(fmtstr,*) '(a',p,',a2,i',l,',a3,a',m,')'
            call fbem_trimall(fmtstr)
            write (selected_unit,fmtstr) trim(errorstr),' (',line_number,'): ',trim(message)
          else
            write(fmtstr,*) '(a',p,',a2,a',n,',a3,a',m,')'
            call fbem_trimall(fmtstr)
            write (selected_unit,fmtstr) trim(errorstr),' (',trim(file_name),'): ',trim(message)
          end if
        end if
      else
        write(fmtstr,*) '(a',p,',a2,a',m,')'
        call fbem_trimall(fmtstr)
        write (selected_unit,fmtstr) trim(errorstr),': ',trim(message)
      end if
      stop
    end if
  end subroutine fbem_error_message

  !! Print a timestamp
  subroutine fbem_timestamp_message(selected_unit,opt)
    implicit none
    ! I/O
    integer           :: selected_unit  !! Unit used to print the message
    integer           :: opt            !! Printing options: 0 (Timestamp="<timestamp>"), 1 ([<timestamp>]\n), 2 (1 following in the same line)
    ! Local
    character(len=8)  :: timestamp_date
    character(len=10) :: timestamp_time
    call date_and_time(timestamp_date,timestamp_time)
    select case (opt)
      case (0)
        write(selected_unit,'(a11,a4,a1,a2,a1,a2,a1,a2,a1,a2,a1,a6,a1)') 'Timestamp="',&
        timestamp_date(1:4),'-',timestamp_date(5:6),'-',timestamp_date(7:8),';',&
        timestamp_time(1:2),':',timestamp_time(3:4),':',timestamp_time(5:10),'"'
      case (1)
        write(selected_unit,'(a1,a4,a1,a2,a1,a2,a1,a2,a1,a2,a1,a6,a1)') '[',&
        timestamp_date(1:4),'-',timestamp_date(5:6),'-',timestamp_date(7:8),';',&
        timestamp_time(1:2),':',timestamp_time(3:4),':',timestamp_time(5:10),']'
      case (2)
        write(selected_unit,'(a1,a4,a1,a2,a1,a2,a1,a2,a1,a2,a1,a6,a1,1x)',advance='no') '[',&
        timestamp_date(1:4),'-',timestamp_date(5:6),'-',timestamp_date(7:8),';',&
        timestamp_time(1:2),':',timestamp_time(3:4),':',timestamp_time(5:10),']'
      case default
        call fbem_error_message(error_unit,0,'fbem_timestamp_message',opt,'opt value not valid')
    end select
  end subroutine fbem_timestamp_message

  !! Print a timestamp with message
  subroutine fbem_timestamp_w_message(selected_unit,opt,message)
    implicit none
    ! I/O
    integer           :: selected_unit  !! Unit used to print the message
    integer           :: opt            !! Printing options: 0 (Timestamp="<timestamp>"), 1 ([<timestamp>]\n), 2 (1 following in the same line)
    character(len=*)  :: message        !! Message to print.
    ! Local
    character(len=fbem_fmtstr) :: fmtstr
    character(len=8)           :: timestamp_date
    character(len=10)          :: timestamp_time
    call date_and_time(timestamp_date,timestamp_time)
    select case (opt)
      case (0)
        write(selected_unit,'(a11,a4,a1,a2,a1,a2,a1,a2,a1,a2,a1,a6,a1)') 'Timestamp="',&
        timestamp_date(1:4),'-',timestamp_date(5:6),'-',timestamp_date(7:8),';',&
        timestamp_time(1:2),':',timestamp_time(3:4),':',timestamp_time(5:10),'"'
      case (1)
        write(selected_unit,'(a1,a4,a1,a2,a1,a2,a1,a2,a1,a2,a1,a6,a1)') '[',&
        timestamp_date(1:4),'-',timestamp_date(5:6),'-',timestamp_date(7:8),';',&
        timestamp_time(1:2),':',timestamp_time(3:4),':',timestamp_time(5:10),']'
      case (2)
        write(selected_unit,'(a1,a4,a1,a2,a1,a2,a1,a2,a1,a2,a1,a6,a1,1x)',advance='no') '[',&
        timestamp_date(1:4),'-',timestamp_date(5:6),'-',timestamp_date(7:8),';',&
        timestamp_time(1:2),':',timestamp_time(3:4),':',timestamp_time(5:10),']'
      case default
        call fbem_error_message(error_unit,0,'fbem_timestamp_message',opt,'opt value not valid')
    end select
    write(fmtstr,*) '(a',len_trim(message),')'
    call fbem_trimall(fmtstr)
    write (selected_unit,fmtstr) trim(message)
  end subroutine fbem_timestamp_w_message

  !! Count the number of words of a string by taking into account that characters between (") are considered a word
  function fbem_count_words(str)
    implicit none
    character(len=*), intent(in)           :: str
    integer                                :: fbem_count_words
    character(len=len(trim(adjustl(str)))) :: s
    integer                                :: i,n
    logical                                :: b,p
    s=trim(adjustl(str))
    n=len_trim(s)
    fbem_count_words=0
    if (n.eq.0) then
      return
    else if (n.eq.1) then
      if (s(1:1).eq.'"') then
        call fbem_error_message(error_unit,0,__FILE__,__LINE__,'incomplete string.')
      else
        fbem_count_words=1
      end if
    else
      b=.true.
      p=.false.
      do i=1,n
        if (i.eq.n) then
          if (p) then
            if (s(i:i).eq.'"') then
            else
              call fbem_error_message(error_unit,0,__FILE__,__LINE__,&
                                      'incomplete string.')
            end if
          else
            if (s(i:i).eq.'"') then
              call fbem_error_message(error_unit,0,__FILE__,__LINE__,&
                                      'incomplete string.')
            else
              if (b) then
                fbem_count_words=fbem_count_words+1
              else
                if (s(i-1:i-1).eq.'"') then
                  call fbem_error_message(error_unit,0,__FILE__,__LINE__,&
                                          'a blank must exist after a string ending.')
                end if
              end if
            end if
          end if
        else
          if (p) then
            if (s(i:i).eq.'"') then
              b=.false.
              p=.false.
            end if
          else
            if (b) then
              if (s(i:i).eq.' ') then
                b=.true.
              else if (s(i:i).eq.'"') then
                fbem_count_words=fbem_count_words+1
                b=.false.
                p=.true.
              else
                fbem_count_words=fbem_count_words+1
                b=.false.
              end if
            else
              if (s(i:i).eq.' ') then
                b=.true.
              else
                if (s(i-1:i-1).eq.'"') then
                  call fbem_error_message(error_unit,0,__FILE__,__LINE__,&
                                          'a blank must exist after a string ending.')
                end if
                b=.false.
              end if
            end if
          end if
        end if
      end do
    end if
  end function fbem_count_words

  !! Select a word in a string by taking into account that characters between (") are considered a word
  function fbem_extract_word(str,sw)
    implicit none
    character(len=*), intent(in)           :: str
    integer                                :: sw
    character(len=len(trim(adjustl(str)))) :: fbem_extract_word
    integer                                :: count_words
    character(len=len(trim(adjustl(str)))) :: s
    integer                                :: i,n
    logical                                :: b,p
    integer                                :: starting,ending
    s=trim(adjustl(str))
    n=len_trim(s)
    count_words=0
    starting=0
    ending=0
    fbem_extract_word=''
    if (n.eq.0) then
      return
    else if (n.eq.1) then
      if (s(1:1).eq.'"') then
        call fbem_error_message(error_unit,0,__FILE__,__LINE__,&
                                'incomplete string.')
      else
        count_words=1
        if (count_words.eq.sw) then
          fbem_extract_word=s(1:1)
          return
        end if
      end if
    else
      b=.true.
      p=.false.
      do i=1,n
        if (i.eq.n) then
          if (p) then
            if (s(i:i).eq.'"') then
              ending=i-1
              if (count_words.eq.sw) then
                fbem_extract_word=adjustl(s(starting:ending))
                return
              end if
            else
              call fbem_error_message(error_unit,0,__FILE__,__LINE__,&
                                      'incomplete string.')
            end if
          else
            if (s(i:i).eq.'"') then
              call fbem_error_message(error_unit,0,__FILE__,__LINE__,&
                                      'incomplete string.')
            else
              if (b) then
                count_words=count_words+1
                if (count_words.eq.sw) then
                  fbem_extract_word=adjustl(s(i:i))
                  return
                end if
              else
                if (s(i-1:i-1).eq.'"') then
                  stop 'Error: a blank must exist after a string ending'
                end if
                ending=i
                if (count_words.eq.sw) then
                  fbem_extract_word=adjustl(s(starting:ending))
                  return
                end if
              end if
            end if
          end if
        else
          if (p) then
            if (s(i:i).eq.'"') then
              ending=i-1
              if (count_words.eq.sw) then
                fbem_extract_word=adjustl(s(starting:ending))
                return
              end if
              b=.false.
              p=.false.
            end if
          else
            if (b) then
              if (s(i:i).eq.' ') then
                b=.true.
              else if (s(i:i).eq.'"') then
                count_words=count_words+1
                starting=i+1
                b=.false.
                p=.true.
              else
                count_words=count_words+1
                starting=i
                b=.false.
              end if
            else
              if (s(i:i).eq.' ') then
                if (s(i-1:i-1).ne.'"') then
                  ending=i-1
                  if (count_words.eq.sw) then
                    fbem_extract_word=s(starting:ending)
                    return
                  end if
                end if
                b=.true.
              else
                if (s(i-1:i-1).eq.'"') then
                  call fbem_error_message(error_unit,0,__FILE__,__LINE__,&
                            'a blank must exist after a string ending.')
                end if
                b=.false.
              end if
            end if
          end if
        end if
      end do
    end if
    if (sw.gt.count_words) then
      write(*,*) trim(str)
      call fbem_error_message(error_unit,0,__FILE__,__LINE__,&
                              'there are less words in the sentence than expected.')
    end if
  end function fbem_extract_word

  !! It finds the first line in a file that matches "[section]", where starting, ending and duplicated blanks are ignored. The file
  !! position remains in the next line.
  subroutine fbem_search_section(selected_unit,section,found)
    implicit none
    ! I/O
    integer                                :: selected_unit !! Unit of the file
    character(len=*)                       :: section       !! Section name to be found. It can't contain starting, ending or duplicated blanks.
    logical                                :: found         !! True if the section has been found, false otherwise.
    ! Local
    character(len=fbem_file_record_length) :: line          ! Line
    integer                                :: word_length   ! Word length
    integer                                :: file_line     ! Current line, at the end, it takes the line where the section line is located.
    integer                                :: ios           ! Error flag
    ! Rewind to the begin of file
    rewind(selected_unit)
    file_line=0
    ! Initialize
    found=.false.
    ! Reading process
    do
      ! Read line
      read(selected_unit,'(a)',iostat=ios) line
      if (is_iostat_end(ios)) then
        exit
      end if
      ! If at the beginning of the file, check if BOM is present.
      if (file_line.eq.0) then
        if ((ichar(line(1:1)).eq.239).and.(ichar(line(2:2)).eq.187).and.(ichar(line(3:3)).eq.191)) line(1:3)=''
      end if
      file_line=file_line+1
      ! Clean the line from spaces at the start and at the end
      call fbem_trim(line)
      ! Check if the line starts with '['
      if (line(1:1).eq.'[') then
        word_length=len_trim(line)
        ! Check if the line ends with ']'
        if (line(word_length:word_length).eq.']') then
          ! Pick up the string between '[' and ']'
          line=line(2:word_length-1)
          ! Clean the string from spaces at the start and at the end and from duplicated blanks
          call fbem_trim2b(line)
          ! Check if the section name is the same
          if (trim(line).eq.trim(section)) then
            found=.true.
            exit
          end if
        end if
      end if
    end do
  end subroutine fbem_search_section

  !! It finds the first line in a file that matches "$section". The file position remains in the next line.
  subroutine fbem_search_section_gmsh(selected_unit,section,found)
    implicit none
    ! I/O
    integer                                :: selected_unit !! Unit of the file
    character(len=*)                       :: section       !! Section name to be found. It can't contain blanks.
    logical                                :: found         !! True if the section has been found, false otherwise.
    ! Local
    character(len=fbem_file_record_length) :: line          ! Line
    integer                                :: word_length   ! Word length
    integer                                :: file_line     ! Current line, at the end, it takes the line where the section line is located.
    integer                                :: ios           ! Error flag
    ! Rewind to the begin of file
    rewind(selected_unit)
    file_line=0
    ! Initialize
    found=.false.
    ! Reading process
    do
      ! Read line
      read(selected_unit,'(a)',iostat=ios) line
      if (is_iostat_end(ios)) then
        exit
      end if
      ! If at the beginning of the file, check if BOM is present.
      if (file_line.eq.0) then
        if ((ichar(line(1:1)).eq.239).and.(ichar(line(2:2)).eq.187).and.(ichar(line(3:3)).eq.191)) line(1:3)=''
      end if
      file_line=file_line+1
      ! Clean the line from spaces at the start and at the end
      call fbem_trim(line)
      ! Check if the line starts with '$'
      if (line(1:1).eq.'$') then
        word_length=len_trim(line)
        ! Pick up the string after '$'
        line=line(2:word_length)
        ! Check if the section name is the same
        if (trim(line).eq.trim(section)) then
          found=.true.
          exit
        end if
      end if
    end do
  end subroutine fbem_search_section_gmsh

  !! It finds a line in a file that starts with "keyword<separator>", where starting, ending and duplicated blanks before the
  !! separator are ignored. The file position remains just after the separator. The finding process stops if a line starting with
  !! "[" is found, i.e. it is done within a section.
  subroutine fbem_search_keyword(selected_unit,keyword,separator,found)
    implicit none
    ! I/O
    integer                                :: selected_unit  !! Unit of the file
    character(len=*)                       :: keyword        !! Keyword to be found. It can't contain starting, ending or duplicated blanks. It can not include the separator.
    character(len=1)                       :: separator      !! Separator between keyword and value, typically '=' or ':'. Only one character.
    logical                                :: found          !! True if the section has been found, false otherwise.
    ! Local
    character(len=fbem_file_record_length) :: line           ! Line
    integer                                :: keyword_length ! Keyword length
    integer                                :: ios            ! Error flag
    integer                                :: i              ! Counter
    logical                                :: doexit         ! Exit control variable
    ! Initialize
    found=.false.
    keyword_length=len_trim(keyword)
    ! Reading process
    do
      read(selected_unit,'(a)',iostat=ios) line
      if (is_iostat_end(ios).eqv.(.true.)) exit
      ! Clean the line from spaces at the start and at the end
      call fbem_trim(line)
      ! Check if the line starts with '[', the start of a section. Then exit.
      if (line(1:1).eq.'[') exit
      ! Erase duplicated blanks
      call fbem_trim2b(line)
      ! Check if it is the keyword
      doexit=.false.
      do i=1,keyword_length
        if (line(i:i).ne.keyword(i:i)) then
          doexit=.true.
          exit
        end if
      end do
      ! If any letter is not the same as the keyword, go to the next line
      if (doexit.eqv.(.true.)) cycle
      ! Check if there is a separator character just after the keyword, if not, go to the next line
      if (line(keyword_length+1:keyword_length+1).ne.separator) then
        if (line(keyword_length+1:keyword_length+1).ne.' ') then
          cycle
        else
          if (line(keyword_length+2:keyword_length+2).ne.separator) cycle
        end if
      end if
      ! The keyword has been found, so the file position must be taken to the previous read.
      found=.true.
      backspace(selected_unit)
      ! Now read until the separator is found
      doexit=.false.
      do
        read(selected_unit,'(a1)',advance='no') line(1:1)
        if (line(1:1).eq.separator) then
          doexit=.true.
          exit
        end if
      end do
      if (doexit.eqv.(.true.)) exit
    end do
  end subroutine fbem_search_keyword

  !! Find a valid unit
  function fbem_get_valid_unit()
    implicit none
    integer :: fbem_get_valid_unit
    integer :: i_unit, min_unit, max_unit
    logical :: exists, already_opened
    fbem_get_valid_unit=output_unit
    min_unit=1
    max_unit=32767
    do i_unit=min_unit,max_unit
      if ((i_unit.ne.5).or.(i_unit.ne.6).or.(i_unit.ne.100).or.(i_unit.ne.101).or.(i_unit.ne.102)) then
        inquire(unit=i_unit,exist=exists,opened=already_opened)
        if ((exists).and.(.not.already_opened)) then
          fbem_get_valid_unit=i_unit
          return
        end if
      end if
    end do
    call fbem_error_message(error_unit,0,'fbem_get_valid_unit',0,'no units available.')
  end function fbem_get_valid_unit

  function fbem_get_current_working_directory()
    use iso_fortran_env
    implicit none
    character(len=4096) :: fbem_get_current_working_directory
    integer :: st
    ! Use GNU Extension, only valid for GNU Fortran
    call getcwd(fbem_get_current_working_directory,st)
    if (st.ne.0) then
      write(*,*) 'getcwd() error with status flag', st
      stop
    end if
    fbem_get_current_working_directory = trim(adjustl(fbem_get_current_working_directory))
  end function fbem_get_current_working_directory

  function fbem_get_os_path_type()
    use iso_fortran_env
    implicit none
    integer :: fbem_get_os_path_type !! 0: unknown, 1: unix-type, 2: Windows-type
    character(len=4096) :: cwd
    cwd = fbem_get_current_working_directory()
    if (cwd(1:1) == '/') then
      fbem_get_os_path_type = 1
    else if ((cwd(2:2) == ':').or.(cwd(1:2).eq.'\\')) then
      fbem_get_os_path_type = 2
    else
      stop 'Unknown OS path type'
    end if
  end function fbem_get_os_path_type

  function fbem_path_is_absolute(filename)
    use iso_fortran_env
    implicit none
    character(len=*), intent(in) :: filename
    logical                      :: fbem_path_is_absolute
    character(len=len(filename)) :: filename_copy
    filename_copy = trim(adjustl(filename))
    fbem_path_is_absolute=.false.
    select case (fbem_get_os_path_type())
      case (1)
        if (scan(filename_copy,'/').eq.1) fbem_path_is_absolute=.true.
      case (2)
        if (len_trim(filename_copy).ge.2) then
          if ((filename_copy(2:2) == ':').or.(filename_copy(1:2).eq.'\\')) then
            fbem_path_is_absolute=.true.
          end if
        end if
    end select
  end function fbem_path_is_absolute

  function fbem_path_is_relative(filename)
    use iso_fortran_env
    implicit none
    character(len=*), intent(in) :: filename
    logical                      :: fbem_path_is_relative
    fbem_path_is_relative = .not.fbem_path_is_absolute(filename)
  end function fbem_path_is_relative

  function fbem_get_dirname(filename)
    use iso_fortran_env
    implicit none
    character(len=*), intent(in) :: filename
    character(len=len(filename)) :: fbem_get_dirname
    character(len=1)             :: separator
    integer                      :: kseparator, kseparator2
    select case (fbem_get_os_path_type())
      case (1)
        kseparator = scan(filename,'/',.true.)
        fbem_get_dirname = ''
        if (kseparator.gt.1) then
           fbem_get_dirname = filename(1:kseparator)
        endif
      case (2)
        kseparator = max(scan(filename,'\',.true.),scan(filename,'/',.true.))
        fbem_get_dirname = ''
        if (kseparator.gt.1) then
           fbem_get_dirname = filename(1:kseparator)
        endif
    end select
  end function fbem_get_dirname

  function fbem_file_exists(filename)
    implicit none
    character(len=*) :: filename
    logical          :: fbem_file_exists
    inquire(file=trim(adjustl(filename)),exist=fbem_file_exists)
  end function fbem_file_exists

! ----------------------------------------------------------------------------------------------------------------------------------
!
! TESTS
!

!  subroutine wordparse_test()
!    implicit none
!    character(len=72) :: line_string
!    character(len=72) :: word
!    integer                               :: line
!    integer                               :: ios
!    integer                               :: nw, kw
!    open(file='test_parser.dat',unit=15,recl=fbem_file_record_length)
!    ios=0
!    line=0
!    do while (ios.eq.0)
!      read(15,'(a)',iostat=ios) line_string
!      line=line+1
!      if (is_iostat_end(ios)) exit
!      nw=fbem_count_words(line_string)
!      write(*,*) 'LINE:',line,'; number of words:',nw
!      do kw=1,nw
!        word=fbem_extract_word(line_string,kw)
!        write(*,*) 'Word ',kw,'=[',trim(word),']'
!      end do
!      write(*,*)
!    end do
!  end subroutine wordparse_test

end module fbem_string_handling
