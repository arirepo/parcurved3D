module mesh_io
  use var_array
  implicit none

  private


  ! high-order mesh data type
  ! consisting arbitrary order 
  ! curved prisms and tets 
  type homesh
     ! private
     ! mpi environment and parallel read/write vars 
     integer :: mpi_rank
     character(len = 800) :: fname

     ! the global (unique) number of each element
     ! n_glob(elem#) = 
     integer, dimension(:), allocatable :: n_glob
     ! the physical coordinates of points in each element
     ! x (:) -> (1:3;xyz, 1:npe) 
     type(double2_array), dimension(:), allocatable :: x
     ! adj(:) -> [adjacent element 1, ..., adjacent element k] 
     type(int_array), dimension(:) , allocatable :: adj
     ! near(:) -> union(node2elem maps(i=1, #vertices)) (global elem #)
     type(int_array), dimension(:) , allocatable :: near

     !                  THE FOLLOWING PART :
     ! ONLY INITIALIZED WHEN NEEDED (NOT APPLICABLE FOR DG)
     !
     ! CONTINUOUS HOMESH FORMAT
     !
     ! npts = the number of final unique interpolation nodes 
     !        (after removing the duplicates nodes in discontinuous grids)
     ! 
     integer :: npts_dg, npts
     real*8 :: tol_dg = 0.0d0
     ! xyz(3, npts_dg) = the physical coordinates of the nodes 
     ! NOTE : npts_dg  is the number of nodes in DG grid (including duplicates)
     !        and we always have: 
     !                              npts_dg > npts
     !     
     real*8, dimension(:, :), allocatable :: xyz
     ! number of boundary (facet) triangles
     integer :: n_bntri 
     ! number of points per boundary triangle 
     integer :: npp_bntri
     ! icon(number of elems, MAX of the number of npe)
     ! bntri_info(1:n_bntri, (npp_bntri+ 3)) = [pt1, pt2, ..., pt_nppbn, tag, 
     !             connected element number , local el. face in that element]
     !
     integer, dimension(:, :), allocatable :: icon, bntri_info
     
   contains

     procedure , public :: init  => init_homesh
     procedure , public :: clean => clean_homesh
     procedure , public :: write => write_homesh
     procedure , public :: read => read_homesh
     ! procedure , public :: append => append_homesh
     procedure, public :: echo => echo_homesh 
     procedure, public :: d2c => convert_discont_2_cont
     procedure, private :: is_pt_in_neigh
     procedure, private :: loc_pt

  end type homesh

  public :: homesh

contains


  ! initializes a high-order mesh object
  !
  subroutine init_homesh(this, mpi_rank, nelem)
    implicit none
    class(homesh), intent(inout) :: this
    integer, intent(in) :: mpi_rank, nelem

    ! mpi environment and parallel read/write vars 
    this%mpi_rank = mpi_rank
    write (this%fname, "(A6,I0.4)") "homesh", this%mpi_rank

    ! the global (unique) number of each element
    allocate(this%n_glob(nelem))
    this%n_glob = 0

    allocate(this%x(nelem))
    allocate(this%adj(nelem))
    allocate(this%near(nelem))

    ! done here
  end subroutine init_homesh

  ! deallocates a homesh object
  subroutine clean_homesh(this)
    implicit none
    class(homesh), intent(inout) :: this

    ! local vars
    integer :: i, nelem

    nelem = size(this%x)

    this%mpi_rank = 0
    if ( allocated(this%n_glob) ) deallocate(this%n_glob)

    do i = 1, nelem
       if ( allocated(this%x(i)%val ) ) deallocate(this%x(i)%val)
       if ( allocated(this%adj(i)%val ) ) deallocate(this%adj(i)%val)
       if ( allocated(this%near(i)%val ) ) deallocate(this%near(i)%val)
    end do
    if ( allocated(this%x) ) deallocate( this%x )
    if ( allocated(this%adj) ) deallocate( this%adj )
    if ( allocated(this%near) ) deallocate( this%near )

    ! 
    this%npts_dg = 0
    this%npts = 0

    if ( allocated(this%xyz) ) deallocate(this%xyz)
    this%n_bntri = 0
    this%npp_bntri = 0
    if ( allocated(this%icon) ) deallocate(this%icon)
    if ( allocated(this%bntri_info) ) deallocate(this%bntri_info)

    ! done here
  end subroutine clean_homesh

  ! binary output for this homesh object
  subroutine write_homesh(this)
    implicit none
    class(homesh), intent(in) :: this

    ! local vars
    integer :: i, unit_id, iostat, nelem

    nelem = size(this%x)

    ! open file ...
    open(newunit=unit_id, file=this%fname, form='unformatted', &
         status='replace', action='write')

    ! write session related vars
    write(unit_id, iostat=iostat) this%mpi_rank
    write(unit_id, iostat=iostat) this%fname

    ! 
    write(unit_id, iostat=iostat) nelem ! number of elements
    write(unit_id, iostat=iostat) this%n_glob
    do i = 1, nelem
       write(unit_id, iostat=iostat) size(this%x(i)%val, 1), size(this%x(i)%val, 2) 
       write(unit_id, iostat=iostat) this%x(i)%val
    end do
    do i = 1, nelem
       write(unit_id, iostat=iostat) size(this%adj(i)%val)
       write(unit_id, iostat=iostat) this%adj(i)%val
    end do
    do i = 1, nelem
       write(unit_id, iostat=iostat) size(this%near(i)%val)
       write(unit_id, iostat=iostat) this%near(i)%val
    end do

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
    integer :: nelem, dim, npe, itmp

    ! bullet proofing
    if ( allocated(this%x) .or. allocated(this%n_glob) ) then
       call this%clean()
    end if

    ! open file ...
    open(newunit=unit_id, file=this%fname, form='unformatted', &
         status='old', action='read')

    ! read session related vars
    read(unit_id, iostat=iostat) this%mpi_rank
    read(unit_id, iostat=iostat) this%fname

    !
    read(unit_id, iostat=iostat) nelem

    allocate( this%n_glob(nelem) )
    read(unit_id, iostat=iostat) this%n_glob

    allocate( this%x(nelem) )
    do i = 1, nelem
       read(unit_id, iostat=iostat) dim, npe
       allocate(this%x(i)%val(dim, npe))
       read(unit_id, iostat=iostat) this%x(i)%val
    end do

    allocate( this%adj(nelem) )
    do i = 1, nelem
       read(unit_id, iostat=iostat) itmp
       allocate(this%adj(i)%val(itmp))
       read(unit_id, iostat=iostat) this%adj(i)%val
    end do

    allocate( this%near(nelem) )
    do i = 1, nelem
       read(unit_id, iostat=iostat) itmp
       allocate(this%near(i)%val(itmp))
       read(unit_id, iostat=iostat) this%near(i)%val
    end do

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
    integer :: i, j, nelem

    print *, '=================================================='
    print *, 'The MPI rank of this chunck of the mesh is : ', this%mpi_rank
    print *, 'The binary file name under which this information' &
         , ' is saved/retrieved is: ', this%fname
    print *, '=================================================='

    ! number of elements
    nelem = size(this%x)

    print *, ' total number of elements in this homesh = ', nelem

    do i = 1, nelem
       print *, 'n_glob(', i, ') = ', this%n_glob(i)
       print *, '---------------------------------------------------------------'
       print *, 'list of points in element # ', i
       print *, '---------------------------------------------------------------'
       do j = 1, size(this%x(i)%val, 2)
          print *, '(', this%x(i)%val(1, j), this%x(i)%val(2, j) &
               , this%x(i)%val(3, j), ')'
       end do
       print *, 'adjacents to elem # ', i , ' are:'
       print *, this%adj(i)%val
       print *, 'elements near elem # ', i , ' are:'
       print *, this%near(i)%val
    end do

    ! done here
  end subroutine echo_homesh

  ! convert a discontinuous format mesh object to a
  ! continuous format usually used in CG FEM solvers
  !
  subroutine convert_discont_2_cont(this, tol_dg)
    implicit none
    class(homesh), intent(inout), target :: this
    real*8, intent(in) :: tol_dg

    ! local vars
    integer :: i, j, tpt, MAX_NPE
    logical :: found_flag
    real*8, pointer :: tx(:, :) => null()

    ! inits
    this%tol_dg = tol_dg

    ! compute the total number of points in the
    ! current discontinuous homesh and max(npes)
    this%npts_dg = 0
    MAX_NPE = size(this%x(1)%val, 2)
    do i = 1, size(this%x) !nelems
       this%npts_dg = this%npts_dg + size(this%x(i)%val, 2)
       if ( MAX_NPE < size(this%x(i)%val, 2) ) then
          MAX_NPE = size(this%x(i)%val, 2)
       end if
    end do

    !
    if ( allocated(this%xyz) ) deallocate(this%xyz)
    allocate( this%xyz(3, this%npts_dg) )

    if ( allocated(this%icon) ) deallocate(this%icon)
    allocate( this%icon( size(this%x) , MAX_NPE) )
    this%icon = -1 ! nice initial value :)

    ! _RESET
    this%npts = 0

    ! find unique points in all elements    
    do i = 1, size(this%x)
       tx => this%x(i)%val
       do j = 1, size(tx, 2) 

          found_flag = .false.
          if ( this%is_pt_in_neigh(xpt = tx(:, j), elnum = i) ) then

             tpt = this%loc_pt(xpt = tx(:, j), x = this%xyz, len = this%npts)

             if ( tpt .ne. -1 ) then ! found; already in there somewhere ...
                found_flag = .true.
             end if

          end if

          if ( .not. found_flag ) then          
             this%npts = this%npts + 1
             this%xyz(:, this%npts) = tx(:, j)
             tpt = this%npts
          end if

          this%icon(i, j) = tpt

       end do

    end do

    ! done here
  end subroutine convert_discont_2_cont

  ! search for a given point <xpt> in a neighborhood
  ! of a given element <elnum> (including all elements 
  ! spatially near the given element) and if found
  ! return .true. otherwise returns .false.
  !
  function is_pt_in_neigh(this, xpt, elnum)
    implicit none
    class(homesh), intent(in), target :: this
    real*8, dimension(:), intent(in) :: xpt
    integer, intent(in) :: elnum
    logical :: is_pt_in_neigh

    ! local vars
    integer :: i, tpt
    integer, pointer :: nears(:) => null()
    real*8, pointer :: tx(:, :) => null()

    nears => this%near(elnum)%val
    is_pt_in_neigh = .false.

    do i = 1, size(nears)
       tx => this%x(nears(i))%val
       tpt = this%loc_pt(xpt = xpt, x = tx , len = size(tx, 2))
       if ( tpt .ne. -1 ) then
          is_pt_in_neigh = .true.
          exit
       end if
    end do

    ! done here
  end function is_pt_in_neigh

  ! locates a point in a given array via distance (L2) norm;
  ! if found within the given tolerance in this%tol_dg
  ! then it returns the location of the point 
  ! in the array; otherwise returns -1
  function loc_pt(this, xpt, x, len)
    implicit none
    class(homesh), intent(in) :: this
    real*8, dimension(:), intent(in) :: xpt
    real*8, dimension(:, :), intent(in) :: x
    integer, intent(in) :: len
    integer :: loc_pt

    ! local vars
    integer :: i
    real*8 :: dx(3), dist

    loc_pt = -1

    do i = 1, len
       dx = x(:, i) - xpt
       dist = sqrt(sum(dx*dx))
       if ( dist > this%tol_dg ) then
          cycle
       else
          loc_pt = i
          return
       end if
    end do

    ! done here
  end function loc_pt

end module mesh_io

! ! a little tester program
! !
! program tester
!   use mesh_io
!   use tetmesher
!   implicit none

!   ! local vars
!   integer :: i, p_pri, d
!   type(homesh) :: thomesh
!   real*8, allocatable :: x(:), y(:), z(:)
!   real*8, allocatable :: x_pri(:, :), x_tet(:, :)

!   ! generate a sample prism
!   p_pri = 8
!   call sample_prism_coords(p = p_pri, shift_x = (/0.0d0, 0.0d0, 0.0d0/) &
!        , x = x_pri)
!   ! generate a sample tet
!   d = p_pri
!   call coord_tet(d = d, x = x , y = y, z = z) 
!   allocate(x_tet(3, size(x)))
!   x_tet(1, :) = x
!   x_tet(2, :) = y
!   x_tet(3, :) = z + 1.0d0

!   ! init high-order mesh object ...
!   call thomesh%init(mpi_rank = 0, nelem = 2)

!   ! add prism(s) 
!   thomesh%n_glob(1) = 1
!   allocate(thomesh%x(1)%val(3, size(x_pri, 2)))
!   thomesh%x(1)%val = x_pri
!   allocate(thomesh%adj(1)%val(5))
!   thomesh%adj(1)%val = (/ -1 , -1, -1, -1, 2 /)
!   allocate(thomesh%near(1)%val(1))
!   thomesh%near(1)%val = (/ 2 /)

!   ! add tet(s)
!   thomesh%n_glob(2) = 2
!   allocate(thomesh%x(2)%val(3, size(x_tet, 2)))
!   thomesh%x(2)%val = x_tet
!   allocate(thomesh%adj(2)%val(4))
!   thomesh%adj(2)%val = (/ 1 , -1, -1, -1 /)
!   allocate(thomesh%near(2)%val(1))
!   thomesh%near(2)%val = (/ 1 /)

!   ! write this mesh
!   call thomesh%write()

!   ! clean it for testing
!   call thomesh%clean()

!   ! read it back again
!   call thomesh%read()

!   ! print on screen
!   call thomesh%echo()

!   ! export to tecplot
!   do i = 1, size(thomesh%x)
!      call export_tet_face_curve(x = thomesh%x(i)%val(1, :) &
!           , y = thomesh%x(i)%val(2, :), z = thomesh%x(i)%val(3, :), mina = 20.0d0 &
!           , maxa = 155.0d0, fname = 'homesh.tec', meshnum = i &
!           , append_it = (i .ne. 1))
!   end do

!   ! convert DG grid to CG grid
!   call thomesh%d2c(tol_dg = 1.0D-14)

!   print *, 'quickly show the CG mesh'
!   print *, '********************************************************'
!   do i = 1, thomesh%npts
!      print *, thomesh%xyz(1, i), thomesh%xyz(2, i), thomesh%xyz(3, i)
!   end do
!   print *, '********************************************************'
 
!   ! finalize
!   call thomesh%clean()

! contains

!   ! generates master elements coords of 
!   ! interpolation points (r,s) and
!   ! stores them in (x,y) 
!   subroutine coord_tri(d, x, y)
!     implicit none
!     integer, intent(in) :: d
!     real*8, dimension(:), allocatable :: x, y

!     ! local vars
!     integer :: npe, i, j, jj
!     real*8 :: dx, dy, xloc, yloc

!     npe = (d+1) * (d+2) / 2
!     dx = 1.0d0 / dble(d)
!     dy = 1.0d0 / dble(d)
!     allocate(x(npe), y(npe))
!     x = 0.0d0; y = 0.0d0

!     jj = 1
!     xloc = 1.0d0 
!     do i = 0, d
!        yloc = 0.0d0
!        do j = 0, i
!           x(jj) = xloc
!           y(jj) = yloc
!           yloc = yloc + dy 
!           jj = jj + 1
!        end do
!        xloc = xloc - dx
!     end do

!     ! done here
!   end subroutine coord_tri

!   subroutine coord_tet(d, x, y, z)
!     implicit none
!     integer, intent(in) :: d
!     real*8, dimension(:), allocatable :: x, y, z

!     ! local vars
!     integer :: npe, i, j, k, jj
!     real*8 :: dx, dy, dz, xloc, yloc, zloc

!     npe = (d+1) * (d+2) * (d+3) / 6
!     dx = 1.0d0 / dble(d)
!     dy = 1.0d0 / dble(d)
!     dz = 1.0d0 / dble(d)

!     allocate(x(npe), y(npe), z(npe))
!     x = 0.0d0; y = 0.0d0; z = 0.0d0
!     jj = 1
!     xloc = 1.0d0 
!     do i = 0, d
!        yloc = 1.0d0 - xloc
!        do j = 0, i
!           zloc = 1.0d0 - xloc - yloc
!           do k = 0, j
!              x(jj) = xloc
!              y(jj) = yloc
!              z(jj) = zloc
!              zloc = zloc - dz
!              jj = jj + 1
!           end do
!           yloc = yloc - dy
!        end do
!        xloc = xloc - dx
!     end do

!     ! done here
!   end subroutine coord_tet

!   ! generates a sample pth-order prism
!   ! for testing mesh_io class
!   subroutine sample_prism_coords(p, shift_x, x)
!     implicit none
!     integer, intent(in) :: p
!     real*8, dimension(:), intent(in) :: shift_x
!     real*8, dimension(:,:), allocatable :: x

!     ! local vars
!     integer :: nz, i, npe, k, jj
!     real*8, dimension(:), allocatable :: x0, y0
!     real*8, dimension((p + 1)) :: z0

!     ! points in the z-directions
!     nz = p + 1

!     ! generate bottom triangle equally spaced point dirstibution
!     call coord_tri(d = p, x = x0, y = y0)
!     z0 = (/ ( (dble(i-1) / dble(nz-1)), i = 1, nz ) /)

!     npe = nz * size(x0)
!     allocate(x(3, npe))

!     jj = 1
!     do i = 1, size(x0)
!        do k = 1, nz
!           x(1, jj) = x0(i) + shift_x(1)
!           x(2, jj) = y0(i) + shift_x(2)
!           x(3, jj) = z0(k) + shift_x(3)
!           jj = jj + 1
!        end do
!     end do

!     ! clean ups
!     if ( allocated(x0) ) deallocate(x0)
!     if ( allocated(y0) ) deallocate(y0)

!     ! done here
!   end subroutine sample_prism_coords

!   ! done here
! end program tester
