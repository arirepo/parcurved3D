module master_elem_distrib
  implicit none

  private


  public :: coord_tri, coord_tet, sample_prism_coords

contains

  ! generates master elements coords of 
  ! interpolation points (r,s) and
  ! stores them in (x,y) 
  subroutine coord_tri(d, x, y)
    implicit none
    integer, intent(in) :: d
    real*8, dimension(:), allocatable :: x, y

    ! local vars
    integer :: npe, i, j, jj
    real*8 :: dx, dy, xloc, yloc

    npe = (d+1) * (d+2) / 2
    dx = 1.0d0 / dble(d)
    dy = 1.0d0 / dble(d)
    allocate(x(npe), y(npe))
    x = 0.0d0; y = 0.0d0

    jj = 1
    xloc = 1.0d0 
    do i = 0, d
       yloc = 0.0d0
       do j = 0, i
          x(jj) = xloc
          y(jj) = yloc
          yloc = yloc + dy 
          jj = jj + 1
       end do
       xloc = xloc - dx
    end do

    ! done here
  end subroutine coord_tri

  subroutine coord_tet(d, x, y, z)
    implicit none
    integer, intent(in) :: d
    real*8, dimension(:), allocatable :: x, y, z

    ! local vars
    integer :: npe, i, j, k, jj
    real*8 :: dx, dy, dz, xloc, yloc, zloc

    npe = (d+1) * (d+2) * (d+3) / 6
    dx = 1.0d0 / dble(d)
    dy = 1.0d0 / dble(d)
    dz = 1.0d0 / dble(d)

    allocate(x(npe), y(npe), z(npe))
    x = 0.0d0; y = 0.0d0; z = 0.0d0
    jj = 1
    xloc = 1.0d0 
    do i = 0, d
       yloc = 1.0d0 - xloc
       do j = 0, i
          zloc = 1.0d0 - xloc - yloc
          do k = 0, j
             x(jj) = xloc
             y(jj) = yloc
             z(jj) = zloc
             zloc = zloc - dz
             jj = jj + 1
          end do
          yloc = yloc - dy
       end do
       xloc = xloc - dx
    end do

    ! done here
  end subroutine coord_tet

  ! generates a sample pth-order prism
  ! for testing mesh_io class
  subroutine sample_prism_coords(p, shift_x, x)
    implicit none
    integer, intent(in) :: p
    real*8, dimension(:), intent(in) :: shift_x
    real*8, dimension(:,:), allocatable :: x

    ! local vars
    integer :: nz, i, npe, k, jj
    real*8, dimension(:), allocatable :: x0, y0
    real*8, dimension((p + 1)) :: z0

    ! points in the z-directions
    nz = p + 1

    ! generate bottom triangle equally spaced point dirstibution
    call coord_tri(d = p, x = x0, y = y0)
    z0 = (/ ( (dble(i-1) / dble(nz-1)), i = 1, nz ) /)

    npe = nz * size(x0)
    allocate(x(3, npe))

    jj = 1
    do i = 1, size(x0)
       do k = 1, nz
          x(1, jj) = x0(i) + shift_x(1)
          x(2, jj) = y0(i) + shift_x(2)
          x(3, jj) = z0(k) + shift_x(3)
          jj = jj + 1
       end do
    end do

    ! clean ups
    if ( allocated(x0) ) deallocate(x0)
    if ( allocated(y0) ) deallocate(y0)

    ! done here
  end subroutine sample_prism_coords

end module master_elem_distrib
