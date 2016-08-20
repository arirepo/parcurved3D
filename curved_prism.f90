module curved_prism
  use var_array
  use tetmesher
  use prism_mesher
  use mpi_comm_mod
  use gen_basis
  use op_cascade
  implicit none

  private

  type brep_interp
     private

     integer :: p ! the order of boundary representation by interpolation
     real*8 :: xv(3, 3), xv_top(3, 3)
     real*8, dimension(:), allocatable :: r, s ! master (r,s) per boundary triangle
     ! the following <x> is the above (r,s) projected on the corresponding CAD face
     ! x(x|y|z, interp. points) 
     real*8, dimension(:, :), allocatable :: x 
     type(basis) :: tbasis ! for Vandermonde related operations

   contains

     procedure, public :: init => init_brep_interp

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
    integer :: p_brep, ttri, tnode
    real*8 :: xv_brep(3, 3), xv_top(3, 3)

    ! mesh outs
    real*8, dimension(:, :), allocatable :: xf
    integer, dimension(:, :), allocatable :: tetcon, neigh
    integer :: nbntri
    integer, dimension(:), allocatable :: bntri
    real*8, allocatable :: uu(:, :)

    ! prism layer vars
    real*8 :: nu = 0.2d0

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
    p_brep = 4
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
       call tbrep(ii)%init(p = p_brep, xv = xv_brep, xv_top = xv_top &
            , tol = tol)
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

    end if

    ! done here
  end subroutine curved_prism_geom

  ! initializes the basis data type
  subroutine init_brep_interp(this, p, xv, xv_top, tol)
    implicit none
    class(brep_interp), intent(inout) :: this
    integer, intent(in) :: p
    real*8, dimension(:,:), intent(in) :: xv, xv_top
    real*8, intent(in) :: tol

    ! the order of this brep
    this%p = p

    ! the coordinates of the bottom boundary facet triangle
    this%xv = xv

    ! the coordinates of the top extruded triangle on the prism
    this%xv_top = xv_top

    ! master point distribution
    call coord_tri(d = this%p, x = this%r, y = this%s)

    ! find the CAD projection 
    call map_master_tri_to_phys(xv = this%xv, r = this%r, s = this%s, tol = tol &
       , x = this%x)

    ! initialize the basis functions
    call this%tbasis%init(this%r, this%s, GEN_TRIANGLE)

    ! done here
  end subroutine init_brep_interp


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

  call curved_prism_geom(tetgen_cmd = 'pq1.414nnY' &
       , facet_file = 'sphere_orient.facet' &
       , cad_file = 'sphere2.iges', nhole = nhole &
       , xh = xh, tol = .03d0, tmpi = tmpi, tvl_info = tvl_info)

  ! nhole = 1
  ! allocate(xh(3))
  ! xh = (/ 10.0d0, 0.0d0, 0.0d0 /)
  ! allocate(tvl_info%tags(40))
  ! tvl_info%tags = (/ (ii, ii = 1, 40) /)

  ! tvl_info%nu = 0.3d0

  ! call curved_prism_geom(tetgen_cmd = 'pq1.214nnY' &
  !      , facet_file = 'civil3_orient.facet' &
  !      , cad_file = 'civil3.iges', nhole = nhole &
  !      , xh = xh, tol = 20.0d0, tmpi = tmpi, tvl_info = tvl_info)



  print *, 'Done! The End!'

  call tmpi%finish()

  ! done here
end program tester
