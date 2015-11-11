module timing
! ---------------------------------------------------
! ///////////////////////////////////////////////////
! This part of this module is written by Taylor Erwin
! 

  implicit none
  private

  integer, parameter :: dp=selected_real_kind(12)

  public :: wtime

  integer  :: wall_cr = 0
  real(dp) :: wall_cr_inv = 0.0_dp

contains

  function wtime()
    integer  :: wtime_int
    real(dp) :: wtime

    continue

    call system_clock(wtime_int)

    if (wall_cr == 0) then
      call system_clock(count_rate=wall_cr)
      wall_cr_inv = 1.0_dp / wall_cr
    end if

    wtime = wtime_int * wall_cr_inv

  end function wtime

! ///////////////////////////////////////////////////
! ---------------------------------------------------

end module timing
