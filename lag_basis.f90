module lag_basis
  ! implements Lagrange basis and their derivatives 
  ! for selected shapes in FEM.
  implicit none

  private


  public :: psi

contains

  subroutine psi(etype, i, r, s, val, der)
    implicit none
    integer, intent(in) :: etype ! element type
    integer, intent(in) :: i     ! node number
    real*8, intent(in)  :: r, s  ! local coords in master elem
    real*8, intent(out) :: val   ! the value of basis  
    real*8, dimension(2), intent(out) :: der   ! the (d/dr,d/ds) of basis

    select case (etype) ! what element?

    case(1) ! triangular elements
       !    2
       !    |\
       !    | \
       !    |  \
       !    3---1
       ! ----------------------------------------------------
       select case (i) ! what node "i" in tri?
       case (1)
          val = 1.0d0-r-s; der(1) = -1.0d0; der(2) = -1.0d0
       case (2)
          val = r; der(1) = 1.0d0; der(2) = 0.0d0
       case (3)
          val = s; der(1) = 0.0d0; der(2) = 1.0d0
       case default
          ! problem
          write(*,*) 'incorrect local node number in tri element! stop.'
          write(*,*) 'error happened in subroutine psi(...).'
          stop
       end select
       ! ----------------------------------------------------
    case(2) ! hex
       ! 4---------3
       ! |         |
       ! |         |
       ! 1---------2
       ! ----------------------------------------------------
       select case (i) ! what node "i" in quad?
       case(1)
          val = 0.25d0 * (1.0d0 - r) * (1.0d0 - s)
          der(1) = -0.25d0 * (1.0d0 - s)
          der(2) = -0.25d0 * (1.0d0 - r)
       case(2)
          val = 0.25d0 * (1.0d0 + r) * (1.0d0 - s)
          der(1) = 0.25d0 * (1.0d0 - s)
          der(2) = -0.25d0 * (1.0d0 + r)
       case(3)
          val = 0.25d0 * (1.0d0 + r) * (1.0d0 + s)
          der(1) = 0.25d0 * (1.0d0 + s)
          der(2) = 0.25d0 * (1.0d0 + r)
       case(4)
          val = 0.25d0 * (1.0d0 - r) * (1.0d0 + s)
          der(1) = -0.25d0 * (1.0d0 + s)
          der(2) =  0.25d0 * (1.0d0 - r)
       case default
          ! problem
          write(*,*) 'incorrect local node number in quad element! stop.'
          write(*,*) 'error happened in subroutine psi(...).'
          stop
       end select
       ! ----------------------------------------------------
    case default
       ! problem : bad element type
       write(*,*) 'incorrect element type in evaluating basis psi(...)! stop.'
       stop
    end select


    ! done here
  end subroutine psi

end module lag_basis
