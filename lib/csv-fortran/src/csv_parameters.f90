!*******************************************************************************
!>
!  Various parameters.

    module csv_parameters

    use csv_kinds

    public

    !! Maximum string length of a real number
    integer(ip),      parameter :: max_real_str_len    = 256
    !! Maximum string length of an integer.
    integer(ip),      parameter :: max_integer_str_len = 256
    !! Default single number format statement
    character(len=*), parameter :: default_sp_fmt   = '(G0)'
    !! Default double number format statement
    character(len=*), parameter :: default_wp_fmt   = '(G0)'
    !! Default quad number format statement
    character(len=*), parameter :: default_qp_fmt   = '(G0)'
    !! Default integer number format statement (for writing integer values to
    !! strings and files).
    character(len=*), parameter :: default_int_fmt  = '(I256)'

    end module csv_parameters
!*******************************************************************************
