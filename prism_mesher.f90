module prism_mesher
  use var_array
  implicit none

  private

  ! information type for viscous layer insertion
  type vl_info
     ! the tags of the boundary facets that need to insert VL
     ! others : will be skipped
     integer, dimension(:), allocatable :: tags

     ! the fraction (min edge length in neighborhood) 
     ! in which boundary facet will be extruded
     real*8 :: nu
  end type vl_info

  public :: comp_bn_tri_normal, find_node2icontag, extrude_bn_tris
  public :: vl_info

contains

  ! computes and stores the normal vector
  ! of each boundary triangle in the direction outward
  ! to the geometry.
  ! This later will be used for extrusion procedure 
  ! and inserting viscous layers.
  !
  ! INPUTS :
  ! npts = number of unrepeated points (nodes)
  ! x(3 * npts) = 3D coords of the points
  ! nquad = num. of quadrilateral surface elements
  ! ntri = num. of triangle surface elems
  ! icontag( (4+1) * nquad + (3+1) * ntri ) = connectivity + tag for each elem
  !  
  ! OUTPUTS:
  ! bn_tri_normal(ntri, 3) = each row has the normal vector of that triangle
  ! 
  subroutine comp_bn_tri_normal(npts, x, nquad, ntri, icontag &
       , bn_tri_normal, min_edg_len)
    implicit none
    integer, intent(in) :: npts
    real*8, dimension(:), intent(in) :: x
    integer, intent(in) :: nquad, ntri
    integer, dimension(:), intent(in), target :: icontag
    real*8, dimension(:, :), allocatable :: bn_tri_normal
    real*8, dimension(:), allocatable :: min_edg_len

    ! local vars
    integer :: i, i1, i2, j, j1, j2
    integer, pointer :: pts(:) => null()
    real*8 :: a(3), b(3), d(3), pts_x(3,3)

    ! bullet proofing ...
    if ( nquad .ne. 0 ) then
       print *, ' the number of boundary quads is not zero! ' &
            , ' Note that comp_bn_tri_normal(...) only works with a ' &
            , ' triangular tesselation of the boundaries. stop'
       stop
    end if

    if ( ntri .eq. 0 ) then
       print *, ' the number of boundary tris is zero! ' &
            , ' fatal error in comp_bn_tri_normal(...). stop'
       stop
    end if

    if ( allocated(bn_tri_normal) ) deallocate(bn_tri_normal)
    allocate(bn_tri_normal(ntri, 3))

    if ( allocated(min_edg_len) ) deallocate(min_edg_len)
    allocate(min_edg_len(ntri))

    ! loop over boundary tris
    do i = 1, ntri
       i1 = (i-1) * 4 + 1
       i2 =  i * 4
       ! get node numbers in this boundary tri
       pts => icontag(i1:(i2-1))
       ! store the coordinates of this boundary triangle
       do j = 1, 3
          j1 = (pts(j)-1) * 3 + 1
          j2 =  pts(j) * 3
          pts_x(j, :) = x(j1:j2)
       end do
       ! form the <a> and <b> vectors consisting
       ! the edges of this triangle ..
       a = pts_x(2,:) - pts_x(1,:)
       b = pts_x(3,:) - pts_x(2,:) 
       d = pts_x(1,:) - pts_x(3,:)

       ! finally compute and store the normal vector
       ! obtained using cross product 
       call comp_normal_cross(a, b, bn_tri_normal(i, :))

       ! store the min edge length for this triangle
       min_edg_len(i) = min( comp_vector_len(a), &
            comp_vector_len(b), comp_vector_len(d) )

    end do ! next boundary triangle ...

    ! done here
  end subroutine comp_bn_tri_normal

  ! computes the vector obtained with
  ! the cross product of <a> and <b> and
  ! stores the result in <c>
  !
  subroutine comp_normal_cross(a, b, c)
    implicit none
    real*8, dimension(:), intent(in) :: a, b
    real*8, dimension(:), intent(out) :: c

    ! performing cross product ...
    c(1) = a(2) * b(3) - b(2) * a(3)
    c(2) = b(1) * a(3) - a(1) * b(3)
    c(3) = a(1) * b(2) - b(1) * a(2)

    ! ! normalize ...
    ! c = c / sqrt(sum(c*c))


    ! done here
  end subroutine comp_normal_cross

  ! computes a vector length
  function comp_vector_len(x)
    implicit none
    real*8, dimension(:), intent(in) :: x
    real*8 :: comp_vector_len

    comp_vector_len = sqrt(sum(x*x))

    ! done here
  end function comp_vector_len

  ! finds the node to triangle connectivity for
  ! the boundary triangles read from a Xpatch facet file
  !
  subroutine find_node2icontag(ntri, icontag, node2icontag)
    implicit none
    integer, intent(in) :: ntri
    integer, dimension(:), intent(in), target :: icontag
    type(int_array), dimension(:) :: node2icontag

    ! local vars
    integer :: i, i1, i2, j
    integer, pointer :: nodes(:) => null()

    do i = 1, ntri
       i1 = 4*(i-1) + 1
       i2 = 4*(i-1) + 3
       nodes => icontag(i1:i2)
       do j = 1, 3
          call push_int_2_array(node2icontag(nodes(j))%val, i)
       end do
    end do

    ! done here
  end subroutine find_node2icontag

  ! extrudes the boundary triangles outward
  ! to a specified extend determined by a 
  ! fraction <nu> of the minimum edge 
  ! length in the neighboring triangles.
  !
  ! NOTE : icontag remains the same and only <x> changes!!!
  ! 
  subroutine extrude_bn_tris(npts, x, bn_tri_normal, icontag, node2icontag &
       ,min_edg_len, tvl_info, x2)
    implicit none
    integer, intent(in) :: npts
    real*8, dimension(:), intent(in) :: x
    real*8, dimension(:,:), intent(in) :: bn_tri_normal
    integer, dimension(:), intent(in) :: icontag
    type(int_array), dimension(:), intent(in) :: node2icontag
    real*8, dimension(:), intent(in) :: min_edg_len
    type(vl_info), intent(in) :: tvl_info
    real*8, dimension(:), allocatable :: x2

    ! local vars
    integer :: i, i1, i2, j, tneigh, tneigh_is_vl, ttag
    real*8 :: a_ave(3), neigh_min, x0(3)

    ! initialize ...
    if ( allocated(x2) ) deallocate(x2)
    allocate(x2(size(x)))
    x2 = x

    ! loop over nodes ...
    do i = 1, npts

       ! locate the coordinates of the current node
       i1 = 3 * (i-1) + 1
       i2 = 3 * i
       x0 = x(i1:i2)

       ! loop over neighboring triangles ...
       a_ave = 0.0d0
       tneigh_is_vl = 0
       do j = 1, size(node2icontag(i)%val)

          ! find the number of the current neighbor triangle
          tneigh = node2icontag(i)%val(j)

          ! find its tag
          ttag = icontag(4*tneigh)

          if ( any(tvl_info%tags .eq. ttag) ) then
             tneigh_is_vl = tneigh_is_vl + 1
          end if
  
          ! find the average of facet normals in the the neighboring triangles
          a_ave = a_ave + bn_tri_normal(tneigh, :)

          ! find the neighborhood minimum edge length 
          if ( j .eq. 1) then
             neigh_min = min_edg_len(tneigh)
          else
             neigh_min = min( neigh_min, min_edg_len(tneigh))
          end if

       end do ! next neighbor

       ! normalize the outward vector 
       a_ave = a_ave / sqrt(sum(a_ave*a_ave))

       ! finally perform the extrusion if is consistent ...
       !
       if ( tneigh_is_vl .eq. 0 ) then 
          ! do nothing; because this point is on a boundary
          ! that does not have viscous layer
       elseif ( tneigh_is_vl .eq. size(node2icontag(i)%val) ) then
          ! all the neighbor facest are on VL region; so
          ! perform extrusion

          ! compute the extrusion
          x2(i1:i2) = x0 + tvl_info%nu * neigh_min * a_ave
       else
          print *, 'The point #', i , ' of given facets is' &
               , ' on the intersection of viscous and inviscid layers! ' &
               , ' These two regions must NOT have an intersection! stop'
          stop
       end if

    end do

    ! done here
  end subroutine extrude_bn_tris

end module prism_mesher
