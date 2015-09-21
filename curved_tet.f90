module curved_tet
  implicit none

  private

  real*8, parameter :: Radius = 2.0d0

  public :: coord_tet, master2curved_tet

contains

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

  subroutine Suv(u, v, S)
    implicit none
    real*8, intent(in) :: u, v
    real*8, dimension(:), intent(out) :: S 

    S(1) = u
    S(2) = v
    S(3) = sqrt(Radius**2 - u**2 - v**2)
    ! S(3) = .5d0 - u - v

    ! done here
  end subroutine Suv

  subroutine master2curved_tet(r,s,t, xA, x, y, z)
    implicit none
    real*8, intent(in) :: r, s, t
    real*8, dimension(:), intent(in) :: xA
    real*8, intent(out) :: x, y, z

    ! local vars
    real*8 :: u, v, alpha
    real*8, dimension(3) :: Sf, xf

    if ( abs(t - 1.0d0) <= 1.0d-15 ) then
       u = r ! simple
       v = s ! simple
    else 
       u = r/(1-t)
       v = s/(1-t)
    end if
    alpha = t

    call Suv(u, v, Sf)

    xf = alpha * xA + (1.0d0 - alpha) * Sf

    x = xf(1)
    y = xf(2)
    z = xf(3)

    ! done here
  end subroutine master2curved_tet

end module curved_tet

program tester
  use tetmesher
  use curved_tet
  implicit none

  ! local vars
  integer :: d, i
  real*8, dimension(:), allocatable :: r, s, t, x, y, z
  real*8, dimension(3) :: xA

  ! tetmeshing vars
  integer :: npts, nquad, ntri, nhole
  real*8, dimension(:), allocatable :: xx, xh 
  integer, dimension(:), allocatable :: icontag
  ! outs
  real*8, dimension(:,:), allocatable :: xf, uu
  integer, dimension(:,:), allocatable :: tetcon, neigh
  integer :: nbntri
  integer, dimension(:), allocatable :: bntri


  ! generate the lagrangian tet. interpolation points
  d = 8
  call coord_tet(d, r, s, t)

  allocate(x(size(r)), y(size(r)), z(size(r))) 

  ! map
  xA = (/ .5d0, .5d0, 3.0d0 /)
  do i = 1, size(r)
     call master2curved_tet(r(i),s(i),t(i), xA, x(i), y(i), z(i))
  end do

  ! show mapped points
  print *, 'x = ', x
  print *, 'y = ', y
  print *, 'z = ', z

  ! generate tetmesh for visualization
  npts = size(r)
  allocate(xx(3 * npts))
  do i = 1, npts
     xx(3*(i-1) + 1) = x(i)
     xx(3*(i-1) + 2) = y(i)
     xx(3*(i-1) + 3) = z(i)
  end do
  nquad = 0
  ntri = 0
  allocate(icontag(0))
  nhole = 0
  allocate(xh(nhole))
  ! 
  call tetmesh('nn', npts, xx, nquad, ntri, icontag, nhole, xh &
       , xf, tetcon, neigh, nbntri, bntri)

  ! write to tecplot
  allocate(uu(1, npts))
  uu = 1.0d0
  call write_u_tecplot_tet(meshnum=1, outfile='cur.tec' &
       , x = xf, icon = tetcon, u = uu, appendit = .false.)

  ! done here
end program tester
