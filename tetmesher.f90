module tetmesher
  use tet_props
  implicit none

  private

    ! local constants (change them for bigger meshes)
    integer, parameter :: max_node = 1000000, max_tet = 1000000, max_triface = 200000
    real*8, parameter :: TOL_SUB_TET = 1.0d-14

  interface
     subroutine tetmesh_intf(command, npts, x, nquad &
          , ntri, icontag, nhole, xh, nptsf, xf, ntet &
          , tetcon, neigh, nbntri, bntri, have_bn_marker) bind( C, name = "main_tetgen_wrapper")
       use iso_c_binding, only : c_int, c_char, c_double
       import
       implicit none
       ! inputs
       character (c_char) :: command(*)
       integer (c_int), value :: npts
       real (c_double) :: x(3 * npts)
       integer (c_int), value :: nquad, ntri 
       integer (c_int) :: icontag(5 * nquad + 4 * ntri)
       integer (c_int), value :: nhole
       real (c_double) :: xh(3 * nhole)

       ! outputs
       integer (c_int) :: nptsf
       real (c_double) :: xf(3 * max_node)
       integer (c_int) :: ntet
       integer (c_int) :: tetcon(4 * max_tet), neigh(4 * max_tet)
       integer(c_int) :: nbntri
       integer(c_int) :: bntri(6 * max_triface)
       integer(c_int), value :: have_bn_marker

     end subroutine tetmesh_intf

  end interface

  public :: read_facet_file, tetmesh
  public :: write_u_tecplot_tet
  public :: export_tet_face_curve
  public :: find_master_elem_sub_tet_conn

contains


  ! read the Xpatch facet file containing 3D surface meshing 
  ! of a complex object in tri/quad mixed format and store 
  ! it in the following 1D arrays:
  !
  ! npts = number of unrepeated points (nodes)
  ! x(3 * npts) = 3D coords of the points
  ! nquad = num. of quadrilateral surface elements
  ! ntri = num. of triangle surface elems
  ! icontag( (4+1) * nquad + (3+1) * ntri ) = connectivity + tag for each elem
  !  
  !
  subroutine read_facet_file(ffile, npts, x, nquad, ntri, icontag)
    implicit none
    character(len = *), intent(in) :: ffile
    integer, intent(out) :: npts
    real*8, dimension(:), allocatable :: x
    integer, intent(out) :: nquad, ntri
    integer, dimension(:), allocatable :: icontag

    ! local vars
    integer :: istat, ii, ntype, npt, tmp1, tmp2
    character(len = 128) :: tline
    integer, dimension(:), allocatable :: itmp1, itmp2

    ! opening the input file
    open ( unit=9, file=ffile , status = 'old', &
         iostat = istat)
    if ( istat /= 0) then 
       print *, 'fatal: could not open XPatch facet <', ffile, '> file! exit'
       stop
    end if

    ! skip the header (four lines)
    do ii = 1, 4
       read(9,*)  
    end do

    ! read number of points
    read(9,*) npts 
    print *, '<', ffile, '> contains ', npts , ' number of points!'
    
    ! store coordinates
    if ( allocated(x) ) deallocate(x)
    allocate( x(3 * npts) )
    do ii = 1, npts
       read(9,*) x(3*ii-2), x(3*ii-1), x(3*ii)
    end do
    print *, 'coordinates are stored from facet file!'

    ! read the types of surface elements
    ! if ntype = 1 then only triangles exist
    ! if ntype = 2 then both quads and triangles (sequentially) exist
    !
    read(9, *) ntype
    select case (ntype)

    case (1) ! then read triangles

       read (9, *) tline
       ! bullet proofing ...
       if (tline .ne. 'Triangles') then
          print *, 'No <Triangles> headline exists in triangle section! stop'
          stop
       end if
       read (9, *) ntri, npt
       print *, '<', ffile, '> contains ', ntri , ' number of triangles!'
       nquad = 0 ! no quad

       if ( allocated(icontag) ) deallocate(icontag)
       npt = npt + 1 ! npt plus tag
       allocate( icontag( npt * ntri ) ) 

       ! fill it
       do ii = 1, ntri
          read(9, *) icontag( ii * npt-3 ), icontag( ii * npt-2 ), icontag( ii * npt-1 ) &
               , tmp1, icontag( ii * npt ), tmp2
       end do

    case (2) ! then read quads and tris sequentially

       read (9, *) tline
       ! bullet proofing ...
       if (tline .ne. 'Quadrilaterals') then
          print *, 'No <Quadrilaterals> headline exists in quad section! stop'
          stop
       end if
       read (9, *) nquad, npt
       print *, '<', ffile, '> contains ', nquad , ' number of quads!'

       if ( allocated(icontag) ) deallocate(icontag)
       npt = npt + 1 ! npt plus tag
       allocate( icontag( npt * nquad ) ) 

       ! fill it
       do ii = 1, nquad
          read(9, *) icontag( ii * npt-4 ), icontag( ii * npt-3 ), icontag( ii * npt-2 ) & 
               , icontag( ii * npt-1 ), tmp1, icontag( ii * npt ), tmp2
       end do

       ! take a temporary copy of <icontag> array
       allocate(itmp1(size(icontag)))
       itmp1 = icontag
       deallocate(icontag)

       ! now proceed to read triangles
       read (9, *) tline
       ! bullet proofing ...
       if (tline .ne. 'Triangles') then
          print *, 'No <Triangles> headline exists in triangle section! stop'
          stop
       end if
       read (9, *) ntri, npt
       print *, '<', ffile, '> contains ', ntri , ' number of triangles!'

       if ( allocated(icontag) ) deallocate(icontag)
       npt = npt + 1 ! npt plus tag
       allocate( icontag( npt * ntri ) ) 

       ! fill it
       do ii = 1, ntri
          read(9, *) icontag( ii * npt-3 ), icontag( ii * npt-2 ), icontag( ii * npt-1 ) &
               , tmp1, icontag( ii * npt ), tmp2
       end do

       ! take a temporary copy of <icontag> array
       allocate(itmp2(size(icontag)))
       itmp2 = icontag
       deallocate(icontag)

       ! concat itmp1 and itmp2 arrays and put the result
       ! in final icontag array to return 
       allocate(icontag(size(itmp1) + size(itmp2)))
       icontag(1:size(itmp1)) = itmp1
       icontag((size(itmp1)+1):(size(itmp1) + size(itmp2))) = itmp2 

    case default
       print *, 'unrecognized types of surface elements in facet file! stop'
       stop

    end select

    ! closing the input file
    close(9)

    ! final clean ups
    if ( allocated(itmp1) ) deallocate(itmp1)
    if ( allocated(itmp2) ) deallocate(itmp2)

    ! done here
  end subroutine read_facet_file

  !
  ! generic tetmesher subroutine
  !
  subroutine tetmesh(cmd, npts, x, nquad, ntri, icontag, nhole, xh &
       , xf, tetcon, neigh, nbntri, bntri, bn_marker)
    use iso_c_binding, only : c_int, c_char, c_double, c_null_char
    implicit none
    character(len = *), intent(in) :: cmd
    integer, intent(in) :: npts, nquad, ntri, nhole
    real*8, dimension(:), intent(in) :: x, xh 
    integer, dimension(:), intent(in) :: icontag
    ! outs
    real*8, dimension(:,:), allocatable :: xf
    integer, dimension(:,:), allocatable :: tetcon, neigh
    integer, intent(out) :: nbntri
    integer, dimension(:), allocatable :: bntri
    integer, intent(in), optional :: bn_marker


    ! local vars
    integer :: ii, leng
    character (kind = c_char, len = 200) :: cmd_in
    integer (c_int) :: npts_in
    real (c_double), allocatable :: x_in(:)
    integer (c_int) :: nquad_in, ntri_in 
    integer (c_int), allocatable :: icontag_in(:)
    integer (c_int) :: nhole_in
    real (c_double), allocatable :: xh_in(:)

    integer (c_int) :: nptsf_out
    real (c_double), allocatable :: xf_out(:)
    integer (c_int) :: ntet_out
    integer (c_int), allocatable :: tetcon_out(:), neigh_out(:)
    integer(c_int) :: nbntri_out
    integer(c_int), allocatable :: bntri_out(:)
    integer :: have_bn_marker
 
    ! assign and init
    cmd_in = cmd//C_NULL_CHAR 
    npts_in = npts
    allocate(x_in(size(x)))
    do ii = 1, size(x)
       x_in(ii) = x(ii)
    end do
    nquad_in = nquad
    ntri_in = ntri
    allocate(icontag_in(size(icontag)))
    do ii = 1, size(icontag)
       icontag_in(ii) = icontag(ii)
    end do
    nhole_in = nhole
    if ( nhole_in > 0 ) then
       allocate(xh_in(size(xh)))
       do ii = 1, size(xh)
          xh_in(ii) = xh(ii)
       end do
    end if

    ! alloc chunck of space for the outputs of C++ function aprriory
    allocate(xf_out(3*max_node), tetcon_out(4*max_tet) &
         , neigh_out(4*max_tet), bntri_out(6*max_triface))

    ! call CPP function
    if ( present(bn_marker) ) then 
       have_bn_marker = bn_marker
    else
       have_bn_marker = 1 ! default
    end if
    call tetmesh_intf(command = cmd_in, npts = npts_in, x = x_in, nquad = nquad_in &
         , ntri = ntri_in, icontag = icontag_in, nhole = nhole_in &
         , xh = xh_in, nptsf = nptsf_out, xf = xf_out, ntet = ntet_out &
         , tetcon = tetcon_out, neigh = neigh_out, nbntri = nbntri_out &
         , bntri = bntri_out, have_bn_marker = have_bn_marker)

    ! set up return values
    !
    ! return nodes
    if ( allocated(xf) ) deallocate(xf)
    allocate(xf(3, nptsf_out))
    xf = reshape( xf_out, (/ 3, nptsf_out /) )

    ! return tets and neighbors
    if ( allocated(tetcon) ) deallocate(tetcon)
    allocate(tetcon(ntet_out, 4))
    tetcon = transpose(reshape( tetcon_out, (/ 4, ntet_out /) ) )

    if ( allocated(neigh) ) deallocate(neigh)
    allocate(neigh(ntet_out, 4))
    neigh = transpose(reshape( neigh_out, (/ 4, ntet_out /) ) )

    ! return boundary triangles
    nbntri = nbntri_out
    leng = 6 * nbntri
    if ( allocated(bntri) ) deallocate(bntri)
    allocate(bntri(leng))
    bntri(1:leng) = bntri_out(1:leng)

    ! clean up
    if ( allocated(x_in) ) deallocate(x_in)
    if ( allocated(icontag_in) ) deallocate(icontag_in)
    if ( allocated(xh_in) ) deallocate(xh_in)
    if ( allocated(xf_out) ) deallocate(xf_out)
    if ( allocated(tetcon_out) ) deallocate(tetcon_out)
    if ( allocated(neigh_out) ) deallocate(neigh_out)
    if ( allocated(bntri_out) ) deallocate(bntri_out)

    ! done here
  end subroutine tetmesh

  ! writes the continious unstructured tetonly grid + solution to 
  ! Tecplot format.
  ! the format of 'u' is assumed to be:  
  !
  ! u(neqs, nnodes)
  !
  subroutine write_u_tecplot_tet(meshnum, outfile, x, icon, u, appendit, is_active)
    implicit none
    integer, intent(in) :: meshnum
    character(len=*), intent(in) :: outfile
    real*8, dimension(:, :), intent(in) :: x
    integer, dimension(:, :), intent(in) :: icon
    real*8, dimension(:,:), intent(in) :: u
    logical, optional :: appendit
    logical, optional :: is_active(:)

    ! local vars
    integer :: i, j, k, neqs, nnodes, nelem

    ! init
    neqs = size(u,1)
    nnodes = size(u,2)
    if ( nnodes .ne. size(x,2) ) then
       print *, 'nnodes .ne. size(x,2)! something is wrong in writing tetmesh! stop'
       stop
    end if

    if ( .not. present(is_active) ) then
       nelem = size(icon, 1)
    else
       nelem = count(is_active)
    end if

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
    write(10, *) 'title = "mesh num. ', meshnum, '"'

    write(10, '(A)', advance = 'no') 'variables = "x", "y", "z"'
    do i = 1, neqs
       write(10, '(A, I1, A)', advance = 'no') ', "u', i,'"'
    end do

    write(10,*) ! new line!


    write(10, '(A, I7, A, I7, A, A)', advance = 'no') 'zone n = ' &
         , nnodes, ', e = ', nelem, ', datapacking = point, ' &
         , 'zonetype = fetetrahedron'
    write(10,*) ! new line!

    ! write coordinates and values of the vector field [u]
    do k = 1, nnodes

       write(10, '(F30.17, A, F30.17, A, F30.17)', advance = 'no') &
            x(1,k), ' ',  x(2,k), ' ', x(3,k)
       do j = 1, neqs
          write(10, '(A, F30.17)', advance = 'no') ' ',  u(j,k)
       end do

       write(10,*)
 
    end do

    ! writing the connectivity matrix
    write(10, *)

    ! write connectivity
    do k = 1, size(icon, 1)
       if ( present(is_active) ) then
          if ( is_active(k) ) then
             write(10, *) ' ',  icon(k,1) &
                  , ' ',  icon(k,2), ' ',  icon(k,3) &
                  , ' ',  icon(k,4)
          else
          end if
       else
          write(10, *) ' ',  icon(k,1) &
               , ' ',  icon(k,2), ' ',  icon(k,3) &
               , ' ',  icon(k,4)
       end if
    end do

    ! close the output file
    close(10)

    ! done here
  end subroutine write_u_tecplot_tet

  ! Performs subtetraheralization of a
  ! general curved tet and filter
  ! tets that have dihedral angle not in range
  ! [mina, maxa] and then finally, write 
  ! the result to tecplot file named fname. 
  ! Has the option to append to the previous export.
  !
  subroutine export_tet_face_curve(x, y, z, mina, maxa, fname &
       , meshnum, append_it, ref_length)
    implicit none
    real*8, dimension(:), intent(in) :: x, y, z
    real*8, intent(in) :: mina, maxa
    character(len = *), intent(in) :: fname
    integer, intent(in) :: meshnum
    logical , intent(in) :: append_it
    real*8, intent(in), optional ::  ref_length

    ! local vars
    integer :: i
    ! tetmeshing vars
    integer :: npts, nquad, ntri, nhole
    real*8, dimension(:), allocatable :: xx, xh 
    integer, dimension(:), allocatable :: icontag
    ! outs
    real*8, dimension(:,:), allocatable :: xf, uu
    integer, dimension(:,:), allocatable :: tetcon, neigh
    integer :: nbntri
    integer, dimension(:), allocatable :: bntri
    ! filter
    logical, dimension(:), allocatable :: is_active

    ! generate tetmesh for visualization
    npts = size(x)
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
    call tetmesh('nQ', npts, xx, nquad, ntri, icontag, nhole, xh &
         , xf, tetcon, neigh, nbntri, bntri, 0)

    ! bullet proof 
    if ( npts .ne. size(xf,2) ) then
       print *, 'Warning : npts .ne. size(xf,2) in export_tet_face_curve(...)!'
       ! print *, 'x = ', x
       ! print *, 'y = ', y
       ! print *, 'z = ', z

       ! print *, 'xf(1, :) = ', xf(1, :)
       ! print *, 'xf(2, :) = ', xf(2, :)
       ! print *, 'xf(3, :) = ', xf(3, :)
       !stop
    end if

    ! filter
    allocate(is_active(size(tetcon, 1)))
    if ( present(ref_length) ) then
       call filter_bad_tets(x = xf, icon = tetcon, mina = mina &
            , maxa= maxa, active = is_active, ref_length = ref_length)
    else
       call filter_bad_tets(x = xf, icon = tetcon, mina = mina &
            , maxa= maxa, active = is_active)
    end if

    ! write to tecplot
    print *, 'writing to Tecplot ...'
    allocate(uu(1, size(xf,2)))
    !
    !
    ! add your desired expression here ...
    if ( all ( y > 0.0d0 ) ) then
       uu = 1.0d0
    else
       uu = 0.0d0
    end if
    !
    !
    !

    call write_u_tecplot_tet(meshnum=meshnum, outfile=fname &
         , x = xf, icon = tetcon, u = uu &
         , appendit = append_it, is_active = is_active)

    print *, 'done writing to Tecplot!'

    ! clean ups
    if ( allocated(xx) ) deallocate(xx)
    if ( allocated(xh) ) deallocate(xh)
    if ( allocated(icontag)) deallocate(icontag)
    if ( allocated(xf) ) deallocate(xf)
    if ( allocated(uu) ) deallocate(uu)
    if ( allocated(tetcon) ) deallocate(tetcon)
    if ( allocated(neigh) ) deallocate(neigh)
    if ( allocated(bntri) ) deallocate(bntri)
    if ( allocated(is_active) ) deallocate(is_active)

    ! done here
  end subroutine export_tet_face_curve

  ! filter bad tetrahedrons
  subroutine filter_bad_tets(x, icon, mina, maxa, active, ref_length)
    implicit none
    real*8, dimension(:, :), intent(in) :: x
    integer, dimension(:, :), intent(in) :: icon
    real*8, intent(in) :: mina, maxa
    logical, dimension(:), intent(out) :: active
    real*8, intent(in), optional :: ref_length

    ! local vars
    integer :: i, j, pt
    real*8 :: xtet(3, 4), ang(6), min_ang, max_ang
    real*8, parameter :: r8_pi = 3.141592653589793D+00

    do i = 1, size(icon, 1) ! loop over tets

       ! fill this tet coords
       do j = 1, 4
          pt = icon(i, j)
          xtet(:, j) = x(:, pt)
       end do

       ! compute measure
       call tetrahedron_dihedral_angles ( tetra = xtet, angle = ang )
       !
       min_ang = minval(abs(ang)) * 180.0d0 / r8_pi
       max_ang = maxval(abs(ang)) * 180.0d0 / r8_pi

       ! decide which one to hide
       if ( (min_ang >= mina) .and. (max_ang <= maxa) ) then
          active(i) = .true.
       else
          active(i) = .false.
       end if

       if ( present( ref_length) ) then
          call tetrahedron_edge_length ( tetra = xtet, edge_length = ang)

          ! decide which one to hide
          if ( (maxval(ang) <= ref_length) ) then
             active(i) = active(i) .and. .true.
          else
             active(i) = .false.
          end if

       end if

    end do

    ! done here
  end subroutine filter_bad_tets

  ! given a master element (r,s,t), the following
  ! subroutine finds the connectivity of sub-tetrahedralization
  ! of that master element. We can use this connectivity
  ! information for visualization purpose
  ! when the master element is transformed to physical
  ! space in the form of a highly curved/twisted element
  ! 
  !
  subroutine find_master_elem_sub_tet_conn(rst, icon)
    implicit none
    real*8, dimension(:, :), intent(in) :: rst
    integer, dimension(:, :), allocatable :: icon

    ! local vars
    integer :: i, npts, i1, i2, nhole, j
    real*8, dimension(:), allocatable :: rst_tmp, xh
    integer, dimension(:), allocatable :: icontag

    ! tetmesher outs
    real*8, dimension(:,:), allocatable :: xf
    integer, dimension(:,:), allocatable :: neigh
    integer :: nbntri
    integer, dimension(:), allocatable :: bntri

    ! permutation
    integer :: ipt, nsub
    integer, dimension(:), allocatable :: iperm
    real*8 :: tetra(3, 4), volume
    integer, dimension(:, :), allocatable :: icon_tmp

    ! size the temporary <rst> buffer
    !
    npts = size(rst, 2)
    allocate(rst_tmp(3 * npts))

    ! now, fill it ...
    do i = 1, npts
       i1 = 3*(i-1) + 1
       i2 = 3*i
       rst_tmp(i1:i2) = rst(:, i)
    end do

    ! nullify the icontag array
    allocate(icontag(0))

    ! nullify holes since we dont have 
    ! any holes in our elements
    nhole = 0
    allocate(xh(nhole))

    ! generate the sub-tetrahedrons
    call tetmesh(cmd = 'nQ', npts = npts, x = rst_tmp &
         , nquad = 0, ntri = 0, icontag = icontag &
         , nhole = nhole, xh = xh &
         , xf = xf, tetcon = icon_tmp &
         , neigh = neigh, nbntri = nbntri, bntri = bntri, bn_marker = 0)

    ! bullet proofing
    if ( size(xf, 2) .ne. size(rst, 2) ) then ! extra points are added then!!!
       print *, 'Steiner points are added in find_master_prism_sub_tet_conn(...)!' &
            , ' stop'
       stop
    end if

    ! match the new <icon> and <xyz> exactly with the input <rst>
    allocate(iperm(size(xf, 2)))

    ! fill the permutation matrix
    do i = 1, size(xf, 2)
       ipt = find_pt( xpt = xf(:, i), x = rst)
       if ( ipt .eq. -1) then
          print *, 'could not find xyz #', i, ' in the original rst!' &
               , ' tetmeshing resulted in missing original points! stop'
          stop
       end if
       iperm(i) = ipt
    end do
    
    ! apply the permutation to correct the connectivity
    call replace_pt(array = icon_tmp, iperm = iperm)

    ! further double check icon_tmp to see zero/negative volume
    ! does not exist at all!
    nsub = 0
    do i = 1, size(icon_tmp, 1)
       do j = 1, 4
          tetra(:, j) = rst(:, icon_tmp(i, j))
       end do
       call tetrahedron_volume ( tetra = tetra, volume = volume)
       if ( volume .le. TOL_SUB_TET) then
          cycle
       else
          nsub = nsub + 1
       end if
    end do

    ! finally, size the final connectivity matrix
    if ( allocated(icon) ) deallocate(icon)
    allocate( icon(nsub, 4) )
    ! and fill it
    nsub = 0
    do i = 1, size(icon_tmp, 1)
       do j = 1, 4
          tetra(:, j) = rst(:, icon_tmp(i, j))
       end do
       call tetrahedron_volume ( tetra = tetra, volume = volume)
       if ( volume .le. TOL_SUB_TET) then
          cycle
       else
          nsub = nsub + 1
          icon(nsub, :) = icon_tmp(i, :)
       end if
    end do


    ! cleanups
    if ( allocated(rst_tmp) ) deallocate(rst_tmp)
    if ( allocated(xh) ) deallocate(xh)
    if ( allocated(icontag) ) deallocate(icontag)
    if ( allocated(xf) ) deallocate(xf)
    if ( allocated(neigh) ) deallocate(neigh)
    if ( allocated(bntri) ) deallocate(bntri)
    if ( allocated(iperm) ) deallocate(iperm)
    if ( allocated(icon_tmp) ) deallocate(icon_tmp)

    ! done here
  end subroutine find_master_elem_sub_tet_conn


  function find_pt(xpt, x)
    implicit none
    real*8, dimension(:), intent(in) :: xpt
    real*8, dimension(:, :), intent(in) :: x
    integer :: find_pt

    ! local vars
    integer :: i

    find_pt = -1

    do i = 1, size(x, 2)
       if (.not. all(x(:, i) .eq. xpt)) then
          cycle
       else
          find_pt = i
          return
       end if
    end do

    ! done here
  end function find_pt

  !
  subroutine replace_pt(array, iperm)
    implicit none
    integer, dimension(:,:), intent(inout) :: array
    integer, dimension(:), intent(in) :: iperm

    ! local vars
    integer :: i, j

    do i = 1, size(array, 1)
       do j = 1, size(array, 2)
             array(i,j) = iperm(array(i,j)) 
       end do
    end do

    ! done here
  end subroutine replace_pt

end module tetmesher

! program tester
!   use tetmesher
!   implicit none


!   ! local vars
!   integer :: ii
!   integer :: npts, nquad, ntri, nhole
!   real*8, dimension(:), allocatable :: x
!   integer, dimension(:), allocatable :: icontag
!   real*8, dimension(:), allocatable :: xh
!   ! outs
!   real*8, dimension(:, :), allocatable :: xf
!   integer, dimension(:, :), allocatable :: tetcon, neigh
!   integer :: nbntri
!   integer, dimension(:), allocatable :: bntri

!   ! Tecplot vars
!   real*8, dimension(:, :), allocatable :: u0

!   ! ! test 1
!   ! call read_facet_file('small_prism.facet', npts, x, nquad, ntri, icontag)
!   ! nhole = 1
!   ! allocate(xh(3))
!   ! xh = (/ 0.11d0, 0.11d0, 0.11d0 /) 

!   ! test 2
!   call read_facet_file('coarse_cylinder.facet', npts, x, nquad, ntri, icontag)
!   nhole = 1
!   allocate(xh(3))
!   xh = (/ 0.75d0, 0.5d0, 0.1d0 /) 

!   ! ! test 3
!   ! call read_facet_file('small_cube.facet', npts, x, nquad, ntri, icontag)
!   ! nhole = 0
!   ! allocate(xh(0))

!   ! !
!   ! ! generic tetmesher subroutine
!   ! !
!   ! subroutine tetmesh(cmd, npts, x, nquad, ntri, icontag, nhole, xh &
!   !      , xf, tetcon, neigh, nbntri, bntri)
!   call tetmesh('pq1.414nn', npts, x, nquad, ntri, icontag, nhole, xh &
!          , xf, tetcon, neigh, nbntri, bntri)

!   ! print outputs
!   print *, 'npts finally generated = ', size(xf, 2)
!   do ii = 1, size(xf, 2)
!      print *, 'node # ', ii, ' x = [', xf(:, ii), ']'
!   end do

!   print *, 'ntet = ', size(tetcon, 1)
!   do ii = 1, size(tetcon, 1)
!      print *, 'tet # ', ii, ' connectivity = [ ', tetcon(ii, :), ']'
!      print *, 'tet # ', ii, ' neighbors = [ ', neigh(ii, :) , ']' 
!   end do
!   print *, 'nbntri = ', nbntri
!   do ii = 1, nbntri
!      print *, 'face tri # ', ii, ' conn = [ ', bntri(6*ii - 5), bntri(6*ii - 4), bntri(6*ii - 3) &
!           , ' ], tag = ', bntri(6*ii - 2), ' adjacent tets = [ ', bntri(6*ii - 1), bntri(6*ii) , ' ]'
!   end do

!   ! export to Tecplot
!   allocate(u0(1, size(xf,2)))
!   u0 = 1.0

!   call write_u_tecplot_tet(meshnum = 1, outfile = 'grd.tec' &
!        , x = xf, icon = tetcon, u = u0, appendit = .false.)

!   ! clean ups
!   if ( allocated(x) ) deallocate(x)
!   if ( allocated(icontag) ) deallocate(icontag)
!   if ( allocated(xh) ) deallocate(xh)
!   if ( allocated(xf) ) deallocate(xf)
!   if ( allocated(tetcon) ) deallocate(tetcon)
!   if ( allocated(neigh) ) deallocate(neigh)
!   if ( allocated(bntri) ) deallocate(bntri)
!   if ( allocated(u0) ) deallocate(u0)

!   ! done here
! end program tester
