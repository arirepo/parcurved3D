module renka_trimesh
  implicit none

  private

  public :: rtrimesh, write_renka_tecplot
  public :: write_renka_tecplot3d

contains


  !
  ! Trimesher for scatter point data in xy plane
  !
  ! "xin(1=x,2=y, :)" is the array of scattered points
  ! "icon" will be reallocated/filled at the end with 
  ! the connectivity matrix of linear triangles.
  ! also "xout' and "yout" will be reallocated/filled.
  !
  subroutine rtrimesh(xin, icon, xout, yout)
    implicit none
    real*8, dimension(:, :), intent(in) :: xin
    integer, dimension(:, :), allocatable :: icon
    real*8, dimension(:), allocatable :: xout, yout

    ! local vars
    integer :: n, ii, ier, ntri, ctris, nadj, pt2, pt3, jj
    integer :: i1, i2, adj_start, adj_num
    integer, dimension(:), allocatable, target :: iadj, iend
    real, dimension(:), allocatable :: x, y, z
    integer, dimension(:), pointer :: adj => null()
    logical, dimension(:), allocatable :: idone
    logical :: bn, add

    ! init
    n = size(xin, 2) ! number of nodes
    allocate(iadj(6*n), iend(n))
    if ( allocated(x)) deallocate(x)
    if ( allocated(y)) deallocate(y)
    if ( allocated(z)) deallocate(z)
    allocate(x(n), y(n), z(n))
    do ii = 1, n
       x(ii) = real(xin(1, ii))
       y(ii) = real(xin(2, ii))
       z(ii) = real(1.0)
    end do

    !
    ! reorder the nodes and data values by sorting on
    !   x-components.  iend is used as temporary storage
    !   for a permutation vector.
    !
    ! call reordr(n, 3, x, y, z, iend)

    !
    ! create the triangulation.  the error flag is ignored
    !   since, after all, nothing can go wrong.
    !
    ier = 0
    call trmesh (n, x, y, iadj, iend, ier)
    if ( ier .ne. 0 ) then
       print *, 'could not make renka triangles! stop'
       stop
    end if

    ! compute number of generated triangles
    ntri = 0
    call compntri (n,6,iadj,iend,0, ntri)
    if ( ntri <= 0 ) then
       print *, 'number of triangles is <= 0 in renka trimesher! stop'
       stop
    end if

    ! init icon
    if ( allocated(icon) ) deallocate(icon)
    allocate( icon(ntri, 3) )
    icon = 0
    ctris = 1 ! current triangles
    allocate( idone(n) )
    idone = .false.

    ! print the adjacency info    
    do ii = 1, n

       ! get the adjancy info for that point ...
       if ( ii .eq. 1) then
          i1 = 1
          i2 = iend(1)
       else
          adj_start = iend(ii-1)+1
          adj_num = iend(ii) - iend(ii-1)
          i1 = adj_start
          i2 = adj_start + adj_num - 1
       end if
       adj => iadj(i1:i2)
       nadj = size(adj)

       ! decide if this node is boundary one
       if ( adj(nadj) .eq. 0 ) then
          bn = .true.
       else
          bn = .false.
       end if

       ! start deciding to add the neigh. triangles to connectivity
       do jj = 1, nadj
          if (jj < nadj) then
             pt2 = adj(jj); pt3 = adj(jj+1)
             if (pt3 .eq. 0) then
                add = .false.
             elseif ( idone(pt2) .or. idone(pt3)) then
                add = .false.
             else
                add = .true.
             end if
          else
             pt2 = adj(nadj); pt3 = adj(1)

             if (pt2 == 0) then
                add = .false.
             elseif ( idone(pt2) .or. idone(pt3) .or. bn ) then
                add = .false.
             else
                add = .true.
             end if
          end if

          if (add) then
             icon(ctris, :) = (/ ii, pt2, pt3 /)
             ctris = ctris + 1
          end if

       end do ! triangles around this point

       idone( ii ) = .true.

    end do ! all points

    ! robust debug
    if ( (ctris-1) .ne. ntri ) then
       print *, 'number of added triangles is not as the analytical number! stop'
       stop
    end if

    ! return values
    if ( allocated( xout) ) deallocate(xout)
    allocate(xout(size(x)))
    xout = x

    if ( allocated( yout) ) deallocate(yout)
    allocate(yout(size(y)))
    yout = y

    ! clean ups
    if ( allocated(iadj) ) deallocate(iadj)
    if ( allocated(iend) ) deallocate(iend)
    if ( allocated(x)) deallocate(x)
    if ( allocated(y)) deallocate(y)
    if ( allocated(z)) deallocate(z)
    if ( associated(adj) ) nullify(adj)
    if ( allocated(idone) ) deallocate(idone)

    ! done here
  end subroutine rtrimesh

  ! writes the unstructured Renka trimesh grid
  ! to Tecplot format 
  !
  subroutine write_renka_tecplot(outfile, icon, x, y, appendit)
    implicit none
    character(len=*), intent(in) :: outfile
    integer, dimension(:, :), intent(in) :: icon
    real*8, dimension(:), intent(in) :: x, y
    logical, intent(in), optional :: appendit

    ! local vars
    integer :: i, j, k

    ! opening for rewrite
    if ( present( appendit ) ) then
       if ( appendit ) then
          open(10, file = outfile, status="old", position="append", action="write")
       else
          open(10, file = outfile, status = 'unknown')
       end if
    else
       open(10, file = outfile, status = 'unknown')
    end if

    ! write header
    write(10, *) 'title = "Renka triangles"'

    write(10, '(A)', advance = 'no') 'variables = "x", "y", "z"'

    write(10,*) ! new line!


    write(10, '(A, I7, A, I7, A, A)', advance = 'no') 'zone n = ' &
         , size(x), ', e = ', size(icon, 1), ', f = fepoint, ' &
         , 'et = triangle'
    write(10,*) ! new line!

    ! write coordinates
    do k = 1, size(x)

       write(10, '(F30.17, A, F30.17)', advance = 'no') &
            x(k), ' ',  y(k)
       write(10, '(A, F30.17)', advance = 'no') ' ',  1.0d0

       write(10,*)
 
    end do

    ! writing the connectivity matrix
    write(10, *)

    do k = 1, size(icon, 1)
       write(10, *) ' ',  icon(k,1) &
            , ' ',  icon(k,2), ' ',  icon(k,3)
    end do

    ! close the output file
    close(10)

    ! done here
  end subroutine write_renka_tecplot

  ! writes the unstructured Renka trimesh 3D surface grid
  ! to Tecplot format 
  ! 
  subroutine write_renka_tecplot3d(outfile, icon, x, y, z, appendit)
    implicit none
    character(len=*), intent(in) :: outfile
    integer, dimension(:, :), intent(in) :: icon
    real*8, dimension(:), intent(in) :: x, y, z
    logical, intent(in), optional :: appendit

    ! local vars
    integer :: i, j, k

    ! opening for rewrite
    if ( present( appendit ) ) then
       if ( appendit ) then
          open(10, file = outfile, status="old", position="append", action="write")
       else
          open(10, file = outfile, status = 'unknown')
       end if
    else
       open(10, file = outfile, status = 'unknown')
    end if

    ! write header
    write(10, *) 'title = "Renka triangles"'

    write(10, '(A)', advance = 'no') 'variables = "x", "y", "z"'

    write(10,*) ! new line!


    write(10, '(A, I7, A, I7, A, A)', advance = 'no') 'zone n = ' &
         , size(x), ', e = ', size(icon, 1), ', f = fepoint, ' &
         , 'et = triangle'
    write(10,*) ! new line!

    ! write coordinates
    do k = 1, size(x)

       write(10, '(F30.17, A, F30.17, A, F30.17)', advance = 'no') &
            x(k), ' ',  y(k), ' ', z(k)

       write(10,*)
 
    end do

    ! writing the connectivity matrix
    write(10, *)

    do k = 1, size(icon, 1)
       write(10, *) ' ',  icon(k,1) &
            , ' ',  icon(k,2), ' ',  icon(k,3)
    end do

    ! close the output file
    close(10)

    ! done here
  end subroutine write_renka_tecplot3d

end module renka_trimesh

! program tester
!   use renka_trimesh
!   implicit none

!   ! local vars
!   integer :: n, n1, n2, ii, jj, jpt
!   real*8, dimension(:,:), allocatable :: x
!   integer, dimension(:, :), allocatable :: icon
!   real*8, dimension(:), allocatable :: xx, yy

!   ! sample scattered node 
!   n1 = 10
!   n2 = 5
!   n = n1 * n2
!   jpt = 1
!   allocate(x(2, n))
!   do ii = 1, n1
!      do jj = 1, n2
!         ! x(1, jpt) = dble(ii-1) / dble(n1-1)
!         ! x(2, jpt) = dble(jj-1) / dble(n2-1)
!         x(1, jpt) = rand()
!         x(2, jpt) = rand()

!         jpt = jpt + 1
!      end do
!   end do

!   ! generate triangulation
!   call rtrimesh(x, icon, xx, yy)

!   ! print the result
!   print *, 'xx = ', xx
!   print *, 'yy = ', yy
!   print *, 'icon = '
!   do n = 1, size(icon, 1)
!      print *, icon(n, :)
!   end do

!   ! write to tecplot file
!   print *, 'writing to the output file grd.dat ...'
!   call write_renka_tecplot(outfile = 'grd.dat', icon = icon &
!        , x = xx, y = yy, appendit = .false.)
!   print *, 'done!'

!   ! clean ups
!   if ( allocated(x) ) deallocate(x)
!   if ( allocated(icon) ) deallocate(icon)
!   if ( allocated(xx) ) deallocate(xx)
!   if ( allocated(yy) ) deallocate(yy)

!   ! done here
! end program tester
