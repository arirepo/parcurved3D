module var_array
  implicit none

  private

  type int_array
     integer, dimension(:), allocatable :: val
  end type int_array


  public :: int_array, push_int_2_array

contains

  subroutine push_int_2_array(a, i)
    implicit none
    integer, dimension(:), allocatable :: a
    integer, intent(in) :: i

    ! local vars
    integer :: n
    integer, dimension(:), allocatable :: itmp

    if ( .not. allocated(a) ) then
       n = 1
    else
       n = size(a) + 1
    end if

    allocate(itmp(n))
    if ( n > 1 ) itmp(1:(n-1)) = a(1:(n-1))
    itmp(n) = i
    call move_alloc(itmp, a)

    ! clean
    if ( allocated(itmp) ) deallocate(itmp)

    ! done here
  end subroutine push_int_2_array

end module var_array
