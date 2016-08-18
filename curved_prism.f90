module curved_prism
  use var_array
  use tetmesher
  use prism_mesher
  use mpi_comm_mod
  implicit none

  private

  public :: curved_prism_geom

contains

  ! 
  subroutine curved_prism_geom(tetgen_cmd, facet_file, cad_file, nhole, xh, tol, tmpi)
    implicit none
    character(len = *), intent(in) :: tetgen_cmd, facet_file, cad_file
    integer, intent(in) :: nhole
    real*8, dimension(:), intent(in) :: xh
    real*8, intent(in) :: tol
    class(mpi_comm_t) :: tmpi

    ! local vars
    !
    !
    ! boundary facet (triangles) vars
    integer :: npts, nquad, ntri
    real*8, dimension(:), allocatable :: x, x2
    integer, dimension(:), allocatable :: icontag
    real*8, dimension(:, :), allocatable :: bn_tri_normal
    real*8, dimension(:), allocatable :: min_edg_len
    type(int_array), dimension(:), allocatable :: node2icontag

    ! mesh outs
    real*8, dimension(:, :), allocatable :: xf
    integer, dimension(:, :), allocatable :: tetcon, neigh
    integer :: nbntri
    integer, dimension(:), allocatable :: bntri
    real*8, allocatable :: uu(:, :)

    ! prism layer vars
    real*8 :: nu = 0.2d0

    if ( tmpi%rank .eq. tmpi%root_rank ) then

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
    call extrude_bn_tris(npts, x, bn_tri_normal, node2icontag,nu,min_edg_len, x2)

    ! !
    ! ! generic tetmesher subroutine
    ! !
    print *, 'generating initial tetmesh of whole domain ...'
    call tetmesh(tetgen_cmd, npts, x2, nquad, ntri, icontag, nhole, xh &
         , xf, tetcon, neigh, nbntri, bntri)
    print *, 'initial tetmesh is done!'

    ! ! export linear tetmesh
    ! allocate(uu(1, size(xf,2)))
    ! uu = 1.0d0
    ! call write_u_tecplot_tet(meshnum=1, outfile='linear_tets.tec', x = xf &
    !      , icon = tetcon, u = uu, appendit = .false.)
    ! if ( allocated(uu) ) deallocate(uu)

    end if

    ! done here
  end subroutine curved_prism_geom


end module curved_prism

program tester
  use curved_prism
  use mpi_comm_mod
  implicit none

  ! local vars
  integer :: nhole
  real*8, allocatable :: xh(:)
  type(mpi_comm_t) :: tmpi

  ! init MPI
  call tmpi%init()


  nhole = 1
  allocate(xh(3))
  xh = 0.0d0

  call curved_prism_geom(tetgen_cmd = 'pq1.414nnY' &
       , facet_file = 'sphere_orient.facet' &
       , cad_file = 'sphere2.iges', nhole = nhole &
       , xh = xh, tol = .03d0, tmpi = tmpi)



  print *, 'Done! The End!'

  call tmpi%finish()

  ! done here
end program tester
