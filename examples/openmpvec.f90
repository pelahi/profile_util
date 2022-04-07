
program OpenMPTest
!---------------------------------------------------------------------
!
!  This is a simple example of using the profiling library
!
!---------------------------------------------------------------------
    use OMP_LIB
    use iso_c_binding, only: C_CHAR, C_NULL_CHAR
    implicit none
    interface 
    ! profile util interfaces 
        subroutine report_binding(string) bind(C)
            use, intrinsic :: iso_c_binding, only: C_CHAR, C_NULL_CHAR
            implicit none
            character(kind=C_CHAR), intent(inout) :: string(*)
        end subroutine
    end interface 
    
    integer, dimension(:), allocatable :: grid
    integer n
    character(len=20000) :: str
    n=1000
    

    call report_binding(str)
    write(*,*) str
    allocate(grid(n))
    deallocate(grid)
end program OpenMPTest

