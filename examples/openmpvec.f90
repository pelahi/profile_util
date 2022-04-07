
program OpenMPTest
!---------------------------------------------------------------------
!
!  This is a simple example of using the profiling library
!
!---------------------------------------------------------------------
    use OMP_LIB
    use iso_c_binding, only: C_CHAR, C_INT, C_NULL_CHAR
    implicit none
    interface 
    ! profile util interfaces 
        subroutine report_parallel_api(string, string_len) bind(C)
            use, intrinsic :: iso_c_binding, only: C_CHAR, C_INT, C_NULL_CHAR
            implicit none
            character(kind=C_CHAR), intent(inout) :: string(*)
            integer(kind=C_INT), intent(inout) :: string_len
        end subroutine
        subroutine report_binding(string, string_len) bind(C)
            use, intrinsic :: iso_c_binding, only: C_CHAR, C_INT, C_NULL_CHAR
            implicit none
            character(kind=C_CHAR), intent(inout) :: string(*)
            integer(kind=C_INT), intent(inout) :: string_len
        end subroutine
    end interface 
    
    integer, dimension(:), allocatable :: grid
    integer :: n, str_len
    character(len=20000) :: str
    n=1000
    

    call report_parallel_api(str, str_len)
    write(*,*) str(str_len)
    call report_binding(str, str_len)
    write(*,*) str(str_len)
    allocate(grid(n))
    deallocate(grid)
end program OpenMPTest

