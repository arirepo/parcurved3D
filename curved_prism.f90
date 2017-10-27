module curved_prism
  use var_array
  use tetmesher
  use prism_mesher
  use mpi_comm_mod
  use gen_basis
  use op_cascade
  use renka_trimesh
  use lag_basis
  use mesh_io
  use master_elem_distrib
  implicit none

  private

  type brep_interp
     private

     integer :: tri_num ! the original number of this bn-tri in the face file
     integer :: p ! the order of boundary representation by interpolation
     real*8 :: xv(3, 3), xv_top(3, 3)
     real*8, dimension(:), allocatable :: r, s ! master (r,s) per boundary triangle
     ! the following <x> is the above (r,s) projected on the corresponding CAD face
     ! x(x|y|z, interp. points) 
     real*8, dimension(:, :), allocatable :: x 

     ! basis functions
     type(basis) :: tbasis ! for Vandermonde related operations
     real*8, dimension(:), allocatable :: psi

     !
     ! boundary triangles visualization and more ...
     integer, dimension(:, :), allocatable :: bn_tris_icon
     real*8, dimension(:, :), allocatable :: x_vis
     logical :: append_flag = .false., enable_vis = .false.
     character(len = 300) :: vis_fname

   contains

     procedure, public :: init => init_brep_interp
     procedure         :: vis_icon => find_bn_tris_icon
     procedure         :: write_bn_tris_tecplot
     procedure         :: rs2xyz_bot => brep_interp_rs2xyz
     procedure, public :: rs2xyz_top => top_interp_rs2xyz
     procedure, public :: gen_prism

  end type brep_interp

  public :: curved_prism_geom

contains

  ! 
  subroutine curved_prism_geom(tetgen_cmd, facet_file, cad_file &
       , nhole, xh, tol, tmpi, tvl_info)
    implicit none
    character(len = *), intent(in) :: tetgen_cmd, facet_file, cad_file
    integer, intent(in) :: nhole
    real*8, dimension(:), intent(in) :: xh
    real*8, intent(in) :: tol
    class(mpi_comm_t) :: tmpi
    type(vl_info), intent(in) :: tvl_info

    ! local vars
    integer :: ii, jj, i1, i2

    !
    ! boundary facet (triangles) vars
    integer :: npts, nquad, ntri
    real*8, dimension(:), allocatable :: x, x2
    integer, dimension(:), allocatable :: icontag
    real*8, dimension(:, :), allocatable :: bn_tri_normal
    real*8, dimension(:), allocatable :: min_edg_len
    type(int_array), dimension(:), allocatable :: node2icontag

    ! boundary representation 
    type(brep_interp), dimension(:), allocatable :: tbrep
    integer, dimension(:), allocatable :: prism_tris
    integer :: ttri, tnode
    real*8 :: xv_brep(3, 3), xv_top(3, 3)

    ! linear mesh outs
    real*8, dimension(:, :), allocatable :: xf
    integer, dimension(:, :), allocatable :: tetcon, neigh
    integer :: nbntri
    integer, dimension(:), allocatable :: bntri
    real*8, allocatable :: uu(:, :)

    ! prism layer vars
    real*8 :: nu = 0.2d0

    ! high-order mesh vars
    integer :: n_homesh
    type(homesh) :: thomesh
    real*8, dimension(:, :), allocatable :: rst_prism
    integer, dimension(:, :), allocatable :: icon_master_prism

    ! viscous layer stretching
    integer :: indx
    real*8, allocatable :: wall(:)

    if ( tmpi%rank .eq. tmpi%root_rank ) then

    ! init CAD file
    call init_IGES_f90(fname = cad_file)

    ! read the facet file
    print *, 'starting curved prisim mesh generator with viscous layer insertation'
    print *, 'reading the facet file ...'
    call read_facet_file(facet_file, npts, x, nquad, ntri, icontag)
    print *, 'the facet file read is complete!'

    print *, 'comp. boundary facet info ...'
    call comp_bn_tri_normal(npts, x, nquad, ntri, icontag &
         , bn_tri_normal, min_edg_len)

    print *, 'find node 2 facet boundary triangles map ...'
    allocate(node2icontag(ntri))
    call find_node2icontag(ntri, icontag, node2icontag)
    print *, 'done with map!'

    print *, 'insert prism layer ...'
    call extrude_bn_tris(npts, x, bn_tri_normal, icontag, node2icontag &
         ,min_edg_len, tvl_info, x2)

    print *, 'find the boundary tris that need to be extruded ...' 
    call detect_prism_tris(ntri, icontag, tvl_info%tags, prism_tris)
    allocate(tbrep(size(prism_tris)))
    print *, 'find the xv_brep and xv_top and initialize the brep object ...'
    do ii = 1, size(prism_tris)
       ttri = prism_tris(ii)
       do jj = 1, 3
          tnode = icontag(4*(ttri-1) + jj)
          i1 = 3 * (tnode-1) + 1
          i2 = 3 * tnode
          xv_brep(:, jj) = x(i1:i2)
          xv_top(:, jj) = x2(i1:i2)
       end do
       print *, ii
       call tbrep(ii)%init(tri_num = ttri, p = tvl_info%p_brep, xv = xv_brep &
            , xv_top = xv_top, tol = tol, enab = tvl_info%enable_bn_tris_vis)
    end do

    ! !
    ! ! generic tetmesher subroutine
    ! !
    print *, 'generating initial tetmesh of whole domain ...'
    call tetmesh(tetgen_cmd, npts, x2, nquad, ntri, icontag, nhole, xh &
         , xf, tetcon, neigh, nbntri, bntri)
    print *, 'initial tetmesh is done!'

    ! export linear tetmesh
    allocate(uu(1, size(xf,2)))
    uu = 1.0d0
    call write_u_tecplot_tet(meshnum=1, outfile='linear_tets.tec', x = xf &
         , icon = tetcon, u = uu, appendit = .false.)
    if ( allocated(uu) ) deallocate(uu)


    ! get the master element point distribution
    ! for the prismatic element of the given order
    !
    call sample_prism_coords(p = 1, shift_x = (/ 0.0d0, 0.0d0, 0.0d0/) &
         , x = rst_prism)
    call find_master_elem_sub_tet_conn(rst = rst_prism, icon = icon_master_prism)

    ! compute the viscous layer stretching
    allocate(wall(tvl_info%nv))
    wall = (/ ((ii - 1.0d0)/(tvl_info%nv-1.0d0) , ii = 1, tvl_info%nv) /)
    wall = wall**2
    indx = 1

    ! init high-order mesh object ...
    n_homesh = size(prism_tris) * (size(wall)-1) ! + number of generated tets
    call thomesh%init(mpi_rank = 0, nelem = n_homesh)


    ! fill the high-order prisms ...
    do ii = 1, size(prism_tris)
       do jj = 1, (size(wall)-1)
          call tbrep(ii)%gen_prism(rst = rst_prism, tb=wall(jj) &
               , tt=wall(jj+1), x = thomesh%x(indx)%val)
          thomesh%n_glob(indx) = indx ! will change when MPI modif. will be added
          allocate(thomesh%adj(indx)%val(1))
          thomesh%adj(indx)%val = indx
          allocate(thomesh%near(indx)%val(1))
          thomesh%near(indx)%val = indx
          indx = indx + 1
       end do
    end do

    ! write this mesh
    call thomesh%write()

    ! export to tecplot
    do ii = 1, size(thomesh%x)

       allocate(uu(1, size(thomesh%x(ii)%val,2)))
       uu = 1.0d0
       call write_u_tecplot_tet(meshnum = ii, outfile = 'homesh2.tec' &
            , x = thomesh%x(ii)%val, icon = icon_master_prism &
            , u = uu, appendit = (ii .ne. 1))
       if ( allocated(uu) ) deallocate(uu)

       ! call export_tet_face_curve(x = thomesh%x(ii)%val(1, :) &
       !      , y = thomesh%x(ii)%val(2, :), z = thomesh%x(ii)%val(3, :) &
       !      , mina = 20.0d0 &
       !      , maxa = 155.0d0, fname = 'homesh.tec', meshnum = ii &
       !      , append_it = (ii .ne. 1))
    end do
    
 
    ! clean it for testing
    call thomesh%clean()


    end if

    ! done here
  end subroutine curved_prism_geom

  ! initializes the basis data type
  subroutine init_brep_interp(this, tri_num, p, xv, xv_top, tol, enab)
    implicit none
    class(brep_interp), intent(inout) :: this
    integer, intent(in) :: tri_num, p
    real*8, dimension(:,:), intent(in) :: xv, xv_top
    real*8, intent(in) :: tol
    logical, intent(in) :: enab

    ! the original number of this bn-tri in the face file
    this%tri_num = tri_num

    ! the order of this brep
    this%p = p

    ! the coordinates of the bottom boundary facet triangle
    this%xv = xv

    ! the coordinates of the top extruded triangle on the prism
    this%xv_top = xv_top

    ! create a name for boundary subtriangluation visualization file
    write (this%vis_fname, "(A8,I0.6)") "bn_tris_", this%tri_num

    ! master point distribution
    call coord_tri(d = this%p, x = this%r, y = this%s)

    ! find the CAD projection 
    call map_master_tri_to_phys(xv = this%xv, r = this%r, s = this%s, tol = tol &
       , x = this%x)

    ! initialize the basis functions
    call this%tbasis%init(this%r, this%s, GEN_TRIANGLE)

    ! reserve an static space for faster
    ! evaluation of the basis functions
    allocate( this%psi(size(this%r)) )
    this%psi = 0.0d0

    ! find the bn-tris connectivity for visualization only
    call this%vis_icon()

    ! write this boundary triangle and its subtriangles
    this%enable_vis = enab
    if ( this%enable_vis ) call this%write_bn_tris_tecplot()

    ! done here
  end subroutine init_brep_interp

  ! finds the connectivity of a triangular mesh
  ! obtained by connecting (this%r, this%s) points
  ! 
  ! NOTE : assuming (r, s) is scattered and has no pattern
  !        a new (r,s) scattering is obtained which leads
  !        to a connectivity. For visualization, these points
  !        must be projected back to CAD face using brep_interp 
  !
  !
  subroutine find_bn_tris_icon(this)
    implicit none
    class(brep_interp), intent(inout) :: this

    ! local vars
    integer :: i
    real*8, dimension(:, :), allocatable :: x_con
    real*8, dimension(:), allocatable :: rtmp, stmp

    allocate(x_con(2, size(this%r)))
    x_con(1, :) = this%r
    x_con(2, :) = this%s 

    call rtrimesh(xin = x_con, icon = this%bn_tris_icon &
         , xout = rtmp, yout = stmp)

    ! project back to interpolated CAD face
    if ( allocated(this%x_vis) ) deallocate(this%x_vis)
    allocate(this%x_vis(3, size(rtmp)))

    do i = 1, size(rtmp)

       call this%rs2xyz_bot(r = rtmp(i), s = stmp(i) &
            , x = this%x_vis(1, i), y = this%x_vis(2, i), z = this%x_vis(3, i))
    end do

    ! clean ups
    if ( allocated(x_con) ) deallocate(x_con)
    if ( allocated(rtmp) ) deallocate(rtmp)
    if ( allocated(stmp) ) deallocate(stmp)

    ! done here
  end subroutine find_bn_tris_icon

  ! writes the sub-triangles of the current 
  ! boundary triangle to tecplot format
  !  
  subroutine write_bn_tris_tecplot(this)
    implicit none
    class(brep_interp), intent(inout) :: this

    !
    call write_renka_tecplot3d(outfile = this%vis_fname, icon = this%bn_tris_icon &
         , x = this%x_vis(1, :), y = this%x_vis(2, :), z = this%x_vis(3, :) &
         , appendit = this%append_flag)

    ! this%append_flag = .true.

    ! done here
  end subroutine write_bn_tris_tecplot

  ! projects a given local (r,s) to the 
  ! physical space using the interpolant
  ! obtained in the current brep_interp object
  !
  subroutine brep_interp_rs2xyz(this, r, s, x, y, z)
    implicit none
    class(brep_interp), intent(inout) :: this
    real*8, intent(in) :: r, s
    real*8, intent(out) :: x, y, z

    ! local vars
    integer :: i
    real*8, dimension(3) :: xyz

    ! evaluate the basis functions 
    call this%tbasis%eval(x0 = r, y0 = s, op = 0, val = this%psi)

    ! compute xyz
    xyz = 0.0d0
    do i = 1, size(this%psi)
       xyz = xyz + this%psi(i) * this%x(:, i)
    end do

    x = xyz(1)
    y = xyz(2)
    z = xyz(3)

    ! done here
  end subroutine brep_interp_rs2xyz

  ! projects a given local (r,s) to the 
  ! physical space using the P1 interpolant
  ! obtained from the top face of the boundary
  ! prismatic frame that contains all sub-prisms
  !
  subroutine top_interp_rs2xyz(this, r, s, x, y, z)
    implicit none
    class(brep_interp), intent(inout) :: this
    real*8, intent(in) :: r, s
    real*8, intent(out) :: x, y, z

    ! local vars
    integer :: i
    real*8, dimension(3) :: xyz
    real*8 :: val   ! the value of basis  
    real*8, dimension(2) :: der   ! the (d/dr,d/ds) of basis

    ! compute xyz
    xyz = 0.0d0
    do i = 1, 3
       call psi(etype = 1, i = i, r = r, s = s, val = val, der = der)
       xyz = xyz + val * this%xv_top(:, i)
    end do

    x = xyz(1)
    y = xyz(2)
    z = xyz(3)

    ! done here
  end subroutine top_interp_rs2xyz

  ! generates a curved boundary conforming prism
  ! 
  subroutine gen_prism(this, rst, tb, tt, x)
    implicit none
    class(brep_interp), intent(inout) :: this
    real*8, dimension(:, :), intent(in) :: rst
    real*8, intent(in) :: tb, tt
    real*8, dimension(:, :), allocatable :: x

    ! local vars
    integer :: i
    real*8 :: r, s, t, xb(3), xt(3)

    if ( allocated(x) ) deallocate(x)
    allocate( x(size(rst,1) , size(rst,2)) )

    do i = 1, size(rst, 2)

       r = rst(1, i)
       s = rst(2, i)
       t = rst(3, i)

       call this%rs2xyz_bot(r = r, s = s, x = xb(1), y = xb(2), z = xb(3))
       call this%rs2xyz_top(r = r, s = s, x = xt(1), y = xt(2), z = xt(3))
       ! map the t value to form a slice
       t = t*tt + (1.0d0 - t) * tb

       x(:, i) = t * xt + (1.0d0 - t) * xb

    end do

    ! done here
  end subroutine gen_prism

end module curved_prism

program tester
  use prism_mesher
  use curved_prism
  use mpi_comm_mod
  implicit none

  ! local vars
  integer :: ii
  integer :: nhole
  real*8, allocatable :: xh(:)
  type(mpi_comm_t) :: tmpi
  type(vl_info) :: tvl_info

  ! init MPI
  call tmpi%init()


  nhole = 1
  allocate(xh(3))
  xh = 0.0d0
  allocate(tvl_info%tags(4))
  tvl_info%tags = (/ 1,2,3,4/)
  tvl_info%nu = 0.3d0
  ! order of boundary representation via polynomials
  tvl_info%p_brep = 4
  tvl_info%enable_bn_tris_vis = .false.
  tvl_info%nv = 10

  call curved_prism_geom(tetgen_cmd = 'pq1.414nnY' &
       , facet_file = 'sphere_orient.facet' &
       , cad_file = 'sphere2.iges', nhole = nhole &
       , xh = xh, tol = .03d0, tmpi = tmpi, tvl_info = tvl_info)

  ! nhole = 1
  ! allocate(xh(3))
  ! xh = (/ 10.0d0, 0.0d0, 0.0d0 /)
  ! allocate(tvl_info%tags(28))
  ! tvl_info%tags = (/ (ii, ii = 1, 28) /)

  ! tvl_info%nu = 0.3d0
  ! tvl_info%p_brep = 1
  ! tvl_info%enable_bn_tris_vis = .false.

  ! call curved_prism_geom(tetgen_cmd = 'pq1.214nnY' &
  !      , facet_file = 'civil3_orient2.facet' &
  !      , cad_file = 'civil3.iges', nhole = nhole &
  !      , xh = xh, tol = 20.0d0, tmpi = tmpi, tvl_info = tvl_info)



  print *, 'Done! The End!'

  call tmpi%finish()

  ! done here
end program tester
