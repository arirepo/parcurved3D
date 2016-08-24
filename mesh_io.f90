module mesh_io
  use var_array
  implicit none

  private


  ! high-order mesh data type
  ! consisting arbitrary order 
  ! curved prisms and tets 
  type homesh

     ! private

     ! number of prisms and number of pts per prism
     integer :: n_pri, npe_pri
     ! the global (unique) number of each prism
     ! n_glob_pri(prism#) = 
     integer, dimension(:), allocatable :: n_glob_pri
     ! the physical coordinates of points in each prism
     ! x_pri (1:3;xyz, 1:npe_pri, 1:n_pri) 
     real*8, dimension(:, :, :), allocatable :: x_pri
     ! icon_pri(prism#, 1:5neighs)
     integer, dimension(:, :) , allocatable :: icon_pri
     ! neigh_pri(1:n_pri) -> union(node2elem maps(i=1, 6vertices))
     type(int_array), dimension(:) , allocatable :: neigh_pri

     ! number of tets and number of pts per tet
     integer :: n_tet, npe_tet
     ! the global (unique) number of each tet
     ! n_glob_tet(tetrahedral#) = 
     integer, dimension(:), allocatable :: n_glob_tet
     ! the physical coordinates of points in each tet
     ! x_tet (1:3;xyz, 1:npe_tet, 1:n_tet) 
     real*8, dimension(:, :, :), allocatable :: x_tet
     ! icon_tet(tetrahedral#, 1:4neighs)
     integer, dimension(:, :) , allocatable :: icon_tet
     ! neigh_tet(1:n_tet) -> union(node2elem maps(i=1, 4vertices))
     type(int_array), dimension(:) , allocatable :: neigh_tet

     ! mpi environment and parallel read/write vars 
     integer :: mpi_rank
     character(len = 800) :: fname


   contains

     procedure , public :: init  => init_homesh
     procedure , public :: clean => clean_homesh
     procedure , public :: write => write_homesh
     procedure , public :: read => read_homesh
     ! procedure , public :: append => append_homesh
     procedure, public :: echo => echo_homesh 

  end type homesh

  public :: homesh

contains


  ! initializes a high-order mesh object
  !
  subroutine init_homesh(this, n_pri, npe_pri, n_tet, npe_tet, mpi_rank)
    implicit none
    class(homesh), intent(inout) :: this
    integer, intent(in) :: n_pri, npe_pri, n_tet, npe_tet, mpi_rank

    ! number of prisms and number of pts per prism
    this%n_pri = n_pri
    this%npe_pri = npe_pri

    ! the global (unique) number of each prism
    ! n_glob_pri(prism#) = 
    allocate(this%n_glob_pri(this%n_pri))
    this%n_glob_pri = 0

    ! the physical coordinates of points in each prism
    ! x_pri (1:3;xyz, 1:npe_pri, 1:n_pri) 
    allocate(this%x_pri(3, this%npe_pri, this%n_pri))
    this%x_pri = 0.0d0

    ! icon_pri(prism#, 1:5neighs)
    allocate(this%icon_pri(this%n_pri, 5))
    this%icon_pri = 0

    ! neigh_pri(1:n_pri) -> union(node2elem maps(i=1, 6vertices))
    allocate(this%neigh_pri(this%n_pri))

    ! number of tets and number of pts per tet
    this%n_tet = n_tet 
    this%npe_tet = npe_tet

    ! the global (unique) number of each tet
    ! n_glob_tet(tetrahedral#) = 
    allocate(this%n_glob_tet(this%n_tet))
    this%n_glob_tet = 0

    ! the physical coordinates of points in each tet
    ! x_tet (1:3;xyz, 1:npe_tet, 1:n_tet) 
    allocate(this%x_tet(3, this%npe_tet, this%n_tet))
    this%x_tet = 0.0d0

    ! icon_tet(tetrahedral#, 1:4neighs)
    allocate(this%icon_tet(this%n_tet, 4))
    this%icon_tet = 0

    ! neigh_tet(1:n_tet) -> union(node2elem maps(i=1, 4vertices))
    allocate(this%neigh_tet(this%n_tet))

    ! mpi environment and parallel read/write vars 
    this%mpi_rank = mpi_rank
    write (this%fname, "(A6,I0.4)") "homesh", this%mpi_rank

    ! done here
  end subroutine init_homesh

  ! deallocates a homesh object
  subroutine clean_homesh(this)
    implicit none
    class(homesh), intent(inout) :: this

    ! local vars
    integer :: i

    ! clean prisms ...
    this%n_pri = 0
    this%npe_pri = 0
    if ( allocated(this%n_glob_pri) ) deallocate(this%n_glob_pri)
    if ( allocated(this%x_pri) ) deallocate(this%x_pri)
    if ( allocated(this%icon_pri) ) deallocate(this%icon_pri)
    do i = 1, size(this%neigh_pri)
       if ( allocated(this%neigh_pri(i)%val) ) deallocate(this%neigh_pri(i)%val)
    end do
    if ( allocated(this%neigh_pri) ) deallocate(this%neigh_pri)

    ! clean tetrahedrons
    this%n_tet = 0
    this%npe_tet = 0
    if ( allocated(this%n_glob_tet) ) deallocate(this%n_glob_tet)
    if ( allocated(this%x_tet) ) deallocate(this%x_tet)
    if ( allocated(this%icon_tet) ) deallocate(this%icon_tet)
    do i = 1, size(this%neigh_tet)
       if ( allocated(this%neigh_tet(i)%val) ) deallocate(this%neigh_tet(i)%val)
    end do
    if ( allocated(this%neigh_tet) ) deallocate(this%neigh_tet)

    ! 
    this%mpi_rank = 0
    ! this%fname = ''

    ! done here
  end subroutine clean_homesh

  ! binary output for this homesh object
  subroutine write_homesh(this)
    implicit none
    class(homesh), intent(in) :: this

    ! local vars
    integer :: i, unit_id, iostat

    ! open file ...
    open(newunit=unit_id, file=this%fname, form='unformatted', &
         status='replace', action='write')

    ! write prisms ...
    write(unit_id, iostat=iostat) this%n_pri, this%npe_pri
    write(unit_id, iostat=iostat) this%n_glob_pri
    write(unit_id, iostat=iostat) this%x_pri
    write(unit_id, iostat=iostat) this%icon_pri
    do i = 1, this%n_pri
       write(unit_id, iostat=iostat) size(this%neigh_pri(i)%val)
       write(unit_id, iostat=iostat) this%neigh_pri(i)%val
    end do

    ! write tets ...
    write(unit_id, iostat=iostat) this%n_tet, this%npe_tet
    write(unit_id, iostat=iostat) this%n_glob_tet
    write(unit_id, iostat=iostat) this%x_tet
    write(unit_id, iostat=iostat) this%icon_tet
    do i = 1, this%n_tet
       write(unit_id, iostat=iostat) size(this%neigh_tet(i)%val)
       write(unit_id, iostat=iostat) this%neigh_tet(i)%val
    end do

    ! write session related vars
    write(unit_id, iostat=iostat) this%mpi_rank
    write(unit_id, iostat=iostat) this%fname

    ! close file ...
    close(unit_id)

    ! done here
  end subroutine write_homesh

  ! binary read for this homesh object
  subroutine read_homesh(this)
    implicit none
    class(homesh), intent(inout) :: this

    ! local vars
    integer :: i, unit_id, iostat
    integer :: itmp

    ! bullet proofing
    if ( allocated(this%x_pri) .or. allocated(this%x_tet) &
         .or. (this%n_pri .ne. 0) .or. (this%n_tet .ne. 0) ) then
       call this%clean()
    end if

    ! open file ...
    open(newunit=unit_id, file=this%fname, form='unformatted', &
         status='old', action='read')

    ! read prisms ...
    read(unit_id, iostat=iostat) this%n_pri, this%npe_pri
    allocate( this%n_glob_pri(this%n_pri) )
    read(unit_id, iostat=iostat) this%n_glob_pri
    allocate( this%x_pri(3, this%npe_pri, this%n_pri) )
    read(unit_id, iostat=iostat) this%x_pri
    allocate( this%icon_pri(this%n_pri, 5) )
    read(unit_id, iostat=iostat) this%icon_pri
    allocate( this%neigh_pri(this%n_pri) )
    do i = 1, this%n_pri
       read(unit_id, iostat=iostat) itmp
       allocate(this%neigh_pri(i)%val(itmp) )
       read(unit_id, iostat=iostat) this%neigh_pri(i)%val
    end do

    ! read tets ...
    read(unit_id, iostat=iostat) this%n_tet, this%npe_tet
    allocate( this%n_glob_tet(this%n_tet) )
    read(unit_id, iostat=iostat) this%n_glob_tet
    allocate( this%x_tet(3, this%npe_tet, this%n_tet) )
    read(unit_id, iostat=iostat) this%x_tet
    allocate( this%icon_tet(this%n_tet, 4) )
    read(unit_id, iostat=iostat) this%icon_tet
    allocate( this%neigh_tet(this%n_tet) )
    do i = 1, this%n_tet
       read(unit_id, iostat=iostat) itmp
       allocate(this%neigh_tet(i)%val(itmp) )
       read(unit_id, iostat=iostat) this%neigh_tet(i)%val
    end do

    ! read session related vars
    read(unit_id, iostat=iostat) this%mpi_rank
    read(unit_id, iostat=iostat) this%fname

    ! close file ...
    close(unit_id)

    ! done here
  end subroutine read_homesh

  ! prints a homesh object on screen
  !
  subroutine echo_homesh(this)
    implicit none
    class(homesh), intent(in) :: this

    ! local vars
    integer :: i, j

    ! print prism(s)
    print *, 'n_pri   = ', this%n_pri
    print *, 'npe_pri = ', this%npe_pri
    do i = 1, this%n_pri
       print *, 'n_glob_pri(', i, ') = ', this%n_glob_pri(i)
    end do
    do i = 1, this%n_pri
       print *, 'In prism # ', i, ' we have:'
       do j = 1, size(this%x_pri, 2)
          print *, 'Point ', j, ' = (', this%x_pri(1, j, i) &
               , this%x_pri(2, j, i), this%x_pri(3, j, i) , ')'
       end do
    end do
    do i = 1, this%n_pri
       print *, 'icon_pri(i, :) = ', this%icon_pri(i, :)
    end do
    do i = 1, this%n_pri
       print *, 'this%neigh_pri(i)%val = ', this%neigh_pri(i)%val 
    end do

    ! print tet(s)
    print *, 'n_tet   = ', this%n_tet
    print *, 'npe_tet = ', this%npe_tet
    do i = 1, this%n_tet
       print *, 'n_glob_tet(', i, ') = ', this%n_glob_tet(i)
    end do
    do i = 1, this%n_tet
       print *, 'In tetrahedron # ', i, ' we have:'
       do j = 1, size(this%x_tet, 2)
          print *, 'Point ', j, ' = (', this%x_tet(1, j, i) &
               , this%x_tet(2, j, i), this%x_tet(3, j, i) , ')'
       end do
    end do
    do i = 1, this%n_tet
       print *, 'icon_tet(i, :) = ', this%icon_tet(i, :)
    end do
    do i = 1, this%n_tet
       print *, 'this%neigh_tet(i)%val = ', this%neigh_tet(i)%val 
    end do

    print *, '=================================================='
    print *, 'The MPI rank of this chunck of the mesh is : ', this%mpi_rank
    print *, 'The binary file name under which this information' &
         , ' is saved/retrieved is: ', this%fname
    print *, '=================================================='

    ! done here
  end subroutine echo_homesh

end module mesh_io

! a little tester program
!
program tester
  use mesh_io
  use tetmesher
  implicit none

  ! local vars
  integer :: p_pri, d
  type(homesh) :: thomesh
  real*8, allocatable :: x(:), y(:), z(:)
  real*8, allocatable :: x_pri(:, :), x_tet(:, :)

  ! generate a sample prism
  p_pri = 1
  call sample_prism_coords(p = p_pri, shift_x = (/0.0d0, 0.0d0, 0.0d0/) &
       , x = x_pri)
  ! generate a sample tet
  d = p_pri
  call coord_tet(d = d, x = x , y = y, z = z) 
  allocate(x_tet(3, size(x)))
  x_tet(1, :) = x
  x_tet(2, :) = y
  x_tet(3, :) = z + 1.0d0

  ! init high-order mesh object ...
  call thomesh%init(n_pri = 1, npe_pri = size(x_pri,2), n_tet = 1 &
       , npe_tet = size(x_tet, 2), mpi_rank = 0)

  ! init prism(s) 
  thomesh%n_glob_pri(1) = 1
  thomesh%x_pri(:, :, 1) = x_pri
  thomesh%icon_pri(1, :) = (/ -1 , -1, -1, -1, 2 /)
  allocate(thomesh%neigh_pri(1)%val(1))
  thomesh%neigh_pri(1)%val(1) = 2

  ! init tet(s)
  thomesh%n_glob_tet(1) = 2
  thomesh%x_tet(:, :, 1) = x_tet
  thomesh%icon_tet(1, :) = (/ 1 , -1, -1, -1/)
  allocate(thomesh%neigh_tet(1)%val(1))
  thomesh%neigh_tet(1)%val(1) = 1


  ! write this mesh
  call thomesh%write()

  ! clean it for testing
  call thomesh%clean()

  ! read it back again
  call thomesh%read()

  ! print on screen
  call thomesh%echo()

  ! export to tecplot
  call export_tet_face_curve(x = thomesh%x_pri(1, :, 1) &
       , y = thomesh%x_pri(2, :, 1), z = thomesh%x_pri(3, :, 1), mina = 20.0d0 &
       , maxa = 155.0d0, fname = 'homesh.tec', meshnum = 1, append_it = .false.)


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

  ! done here
end program tester
