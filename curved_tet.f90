module curved_tet
  use tetmesher
  use tet_props
  use lag_basis
  use op_cascade
  implicit none

  private

  real*8, parameter :: Radius = 2.0d0

  ! type int_array
  !    integer, dimension(:), allocatable :: val
  ! end type int_array

  public :: coord_tet, master2curved_tet
  public :: export_tet_face_curve, master2curved_edg_tet
  public :: curved_tetgen_geom

  ! testers
  public :: tester1

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

  subroutine master2curved_tet(r,s,t, uv, xA, x, y, z)
    implicit none
    real*8, intent(in) :: r, s, t
    real*8, dimension(:, :), intent(in) :: uv
    real*8, dimension(:), intent(in) :: xA
    real*8, intent(out) :: x, y, z

    ! local vars
    integer :: ii
    real*8 :: u, v, alpha
    real*8, dimension(3) :: Sf, xf
    !
    real*8 :: val   ! the value of basis  
    real*8, dimension(2) :: der, uv_fin 


    if ( abs(t - 1.0d0) <= 1.0d-15 ) then
       u = r ! simple
       v = s ! simple
    else 
       u = r/(1-t)
       v = s/(1-t)
    end if
    alpha = t

    ! compute final uv
    uv_fin = 0.0d0
    do ii = 1, 3
       call psi(etype = 1, i = ii, r = u, s = v, val = val, der = der)
       uv_fin = uv_fin + val * uv(:, ii)
    end do

    ! compute surface points
    call Suv(u = uv_fin(1), v= uv_fin(2), S = Sf)

    xf = alpha * xA + (1.0d0 - alpha) * Sf

    x = xf(1)
    y = xf(2)
    z = xf(3)

    ! done here
  end subroutine master2curved_tet

  subroutine master2curved_edg_tet(r,s,t, uv, x3, x4, x, y, z)
    implicit none
    real*8, intent(in) :: r, s, t
    real*8, dimension(:, :), intent(in) :: uv ! uv(1:2, 1:2)
    real*8, dimension(:), intent(in) :: x3, x4
    real*8, intent(out) :: x, y, z

    ! local vars
    real*8 :: u, v, alpha, u0
    real*8, dimension(3) :: Sf, xf, Sc
    !
    real*8, dimension(2) :: uv_fin 


    if ( abs(t - 1.0d0) <= 1.0d-15 ) then
       u = r ! simple
       v = s ! simple
    else 
       u = r/(1-t)
       v = s/(1-t)
    end if
    alpha = t

    ! compute final uv
    if ( abs(v - 1.0d0) <= 1.0d-15 ) then
       u0 = u ! simple
    else 
       u0 = u/(1-v)
    end if

    uv_fin = uv(:, 1) + u0 * (uv(:, 2) - uv(:, 1)) 

    ! compute surface points
    call Suv(u = uv_fin(1), v= uv_fin(2), S = Sc)

    Sf = v * x3 + (1.0d0 - v) * Sc

    xf = alpha * x4 + (1.0d0 - alpha) * Sf

    x = xf(1)
    y = xf(2)
    z = xf(3)

    ! done here
  end subroutine master2curved_edg_tet

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
    call tetmesh('nn', npts, xx, nquad, ntri, icontag, nhole, xh &
         , xf, tetcon, neigh, nbntri, bntri, 0)

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
    allocate(uu(1, npts))
    uu = 1.0d0
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

  !
  ! This tests a set of three tetrahedrons
  ! that exhibit two possible cases of mounting
  ! on a curved surface:
  !
  ! 1 - one face on surface
  ! 2 - only one edge on the surface
  ! 
  ! a sphere is used as exact surface representation
  ! 
  subroutine tester1()
    implicit none

    ! local vars
    integer :: d, i
    real*8, dimension(:), allocatable :: r, s, t, x, y, z
    real*8, dimension(3) :: xA, x3, x4
    real*8, dimension(2, 3) :: uv

    ! generate the lagrangian tet. interpolation points
    d = 8
    call coord_tet(d, r, s, t)
    allocate(x(size(r)), y(size(r)), z(size(r))) 

    ! element 1
    xA = (/ .15d0, .15d0, 2.2d0 /)
    uv = reshape( (/ 0.0d0, 0.0d0, .5d0, 0.0d0 &
         , 0.0d0, 0.5d0 /), (/2, 3/) )

    do i = 1, size(r)
       call master2curved_tet( r = r(i), s = s(i), t = t(i), uv = uv &
            , xA = xA, x = x(i), y = y(i), z = z(i) )

       ! call master2curved_tet(r(i),s(i),t(i), xA, x(i), y(i), z(i))
    end do

    ! export the generated curved element
    call export_tet_face_curve(x = x, y=y, z=z, mina = 20.0d0 &
         , maxa = 155.0d0, fname = 'curved.tec', meshnum = 1, append_it = .false.)  

    ! element 2
    xA = (/ .38d0, .38d0, 2.2d0 /)
    uv = reshape( (/ 0.5d0, 0.0d0, 0.5d0, 0.5d0 &
         , 0.0d0, 0.5d0 /), (/2, 3/) )

    do i = 1, size(r)
       call master2curved_tet( r = r(i), s = s(i), t = t(i), uv = uv &
            , xA = xA, x = x(i), y = y(i), z = z(i) )

       ! call master2curved_tet(r(i),s(i),t(i), xA, x(i), y(i), z(i))
    end do

    ! export the generated curved element
    call export_tet_face_curve(x = x, y=y, z=z, mina = 20.0d0 &
         , maxa = 155.0d0, fname = 'curved.tec', meshnum = 2, append_it = .true.)  

    ! element 3
    xA = (/ .22d0, .22d0, 3.0d0 /)
    uv = reshape( (/ 0.5d0, 0.0d0, 0.0d0, 0.5d0 &
         , 0.0d0, 0.0d0 /), (/2, 3/) )
    x3 = (/ .15d0, .15d0, 2.2d0 /)
    x4 = (/ .38d0, .38d0, 2.2d0 /)

    do i = 1, size(r)
       call master2curved_edg_tet( r = r(i), s = s(i), t = t(i), uv = uv &
            , x3 = x3, x4 = x4, x = x(i), y = y(i), z = z(i) )
    end do

    ! export the generated curved element
    call export_tet_face_curve(x = x, y=y, z=z, mina = 20.0d0 &
         , maxa = 155.0d0, fname = 'curved.tec', meshnum = 3, append_it = .true.)  

    ! clean ups
    if ( allocated(r) ) deallocate(r)
    if ( allocated(s) ) deallocate(s)
    if ( allocated(t) ) deallocate(t)
    if ( allocated(x) ) deallocate(x)
    if ( allocated(y) ) deallocate(y)
    if ( allocated(z) ) deallocate(z)

    ! done here
  end subroutine tester1

  ! 
  subroutine curved_tetgen_geom(tetgen_cmd, facet_file, cad_file, nhole, xh, tol)
    implicit none
    character(len = *), intent(in) :: tetgen_cmd, facet_file, cad_file
    integer, intent(in) :: nhole
    real*8, dimension(:), intent(in) :: xh
    real*8, intent(in) :: tol

    ! local vars
    integer :: ii, jj
    integer :: npts, nquad, ntri
    real*8, dimension(:), allocatable :: x
    integer, dimension(:), allocatable :: icontag

    ! outs
    real*8, dimension(:, :), allocatable :: xf
    integer, dimension(:, :), allocatable :: tetcon, neigh
    integer :: nbntri
    integer, dimension(:), allocatable :: bntri
    ! integer, dimension(:, :), allocatable :: bntri2bntri
    ! type(int_array), dimension(:), allocatable :: node2bntri

    real*8, allocatable :: uu(:, :)

    ! domain decomposition (graph partitioning) vars
    integer, dimension(:), allocatable :: xadj, adj, vwgt, part
    integer :: nparts
    logical, dimension(:), allocatable :: vis_mask

    ! CAD corresponding data struct
    integer, allocatable :: cent_cad_found(:) !nbntri
    real*8, allocatable :: uvc(:)
    integer, dimension(:, :), allocatable :: tet_shifted
    integer, dimension(:), allocatable :: tet2bn_tri
    ! integer :: neigh_CAD_face(3)
    ! logical :: is_CAD_bn_tri

    ! master element data struct
    integer :: dd, indx, CAD_face
    real*8, dimension(:), allocatable :: rr, ss, tt, xx, yy, zz
    real*8, dimension(3) :: xA, x3, x4
    real*8, dimension(2, 3) :: uv
    real*8, dimension(3, 3) :: xbot

    ! visualization data struct
    real*8 :: xtet(3, 4), lens(6)
    real*8 :: ref_length

    ! read the facet file
    print *, 'starting curved tetrahedral mesh generator'
    print *, 'reading the facet file ...'
    call read_facet_file(facet_file, npts, x, nquad, ntri, icontag)
    print *, 'the facet file read is complete!'

    ! !
    ! ! generic tetmesher subroutine
    ! !
    print *, 'generating initial tetmesh of whole domain ...'
    call tetmesh(tetgen_cmd, npts, x, nquad, ntri, icontag, nhole, xh &
         , xf, tetcon, neigh, nbntri, bntri)
    print *, 'initial tetmesh is done!'

    ! ! find the boundary tri connectivity map
    ! ! useful for speedup the code when deciding on
    ! ! UV-projection or closest point
    ! call find_bntri2bntri_map(nbntri = nbntri, bntri = bntri, bntri2bntri = bntri2bntri)

    ! ! bullet proofing ...
    ! if ( any ( bntri2bntri .eq. -1) ) then
    !    print *, 'boundary triangles are not all connected together! stop'
    !    stop
    ! end if

    ! allocate(node2bntri(size(xf, 2)))
    ! call find_node2bntri_map(nbntri = nbntri, bntri = bntri, node2bntri = node2bntri)


    ! export linear tetmesh
    allocate(uu(1, size(xf,2)))
    uu = 1.0d0
    call write_u_tecplot_tet(meshnum=1, outfile='linear_tets.tec', x = xf &
         , icon = tetcon, u = uu, appendit = .false.)
    if ( allocated(uu) ) deallocate(uu)

    call find_tet_adj_csr_zero_based(neigh = neigh, xadj = xadj, adj = adj)

    nparts = 10
    allocate(vwgt(size(neigh, 1)), part(size(neigh, 1)))
    vwgt = 1
    call call_metis_graph_parti(xadj = xadj, adjncy = adj, vwgt = vwgt &
         , nparts = nparts, part = part)

    ! write to tecplot
    allocate(vis_mask(size(neigh, 1)))
    allocate(uu(1, size(xf,2)))
    do ii = 1, nparts

       uu = dble(ii)
       vis_mask = (part .eq. (ii-1))

       if ( ii .eq. 1) then
          call write_u_tecplot_tet(meshnum=ii, outfile='partitioned.tec', x = xf &
               , icon = tetcon, u = uu, appendit = .false., is_active = vis_mask)
       else
          call write_u_tecplot_tet(meshnum=ii, outfile='partitioned.tec', x = xf &
               , icon = tetcon, u = uu, appendit = .true., is_active = vis_mask)
       end if

    end do
    if ( allocated(uu) ) deallocate(uu)
    if ( allocated(vis_mask) ) deallocate(vis_mask)

    ! find the CAD face of boundary triangles
    call find_bn_tris_CAD_face(cad_file = cad_file, nbntri = nbntri &
         , bntri = bntri, xf = xf, cent_cad_found = cent_cad_found &
         , uvc = uvc, tol = tol)

    ! shift tetcon for each tet such that the first three nodes
    ! are matching the boundary triangle face. This is required
    ! before we apply our analytical transformation
    !
    allocate( tet_shifted(size(tetcon, 1), size(tetcon, 2)) )
    allocate( tet2bn_tri(size(tetcon, 1)) )
    call shift_tetcon(nbntri = nbntri, bntri = bntri &
         , tetcon = tetcon, tetcon2 = tet_shifted, tet2bn_tri = tet2bn_tri)
 
    ! print *, 'tet2bn_tri = ', tet2bn_tri
    ! /////////
    ! /////////
    ! generate the lagrangian tet. interpolation points
    dd = 6
    indx = 1
    call coord_tet(dd, rr, ss, tt)
    allocate(xx(size(rr)), yy(size(rr)), zz(size(rr))) 

    ! map one face curved tets
    do ii = 1, size(tet2bn_tri)
       if ( tet2bn_tri(ii) .eq. -1 ) cycle !interior tet
       ! 
       CAD_face = cent_cad_found(tet2bn_tri(ii))
       ! print *, 'CAD_face =' , CAD_face
       if ( CAD_face .eq. -1 ) cycle !bn tet not on CAD database

       !
       xA = xf(:, tet_shifted(ii, 4))

       ! compute local uv on that face
       do jj = 1, 3
          call xyz2uv_f90(CAD_face = CAD_face &
               , xyz = xf(:, tet_shifted(ii, jj)), uv = uv(:,jj), tol = 1.0d-6)
       end do

       ! compute xbot
       do jj = 1, 3
          xbot(:, jj) = xf(:, tet_shifted(ii, jj))
       end do

       ! ! compute the CAD face tag of the neighbors to this boundary
       ! ! triangle. 
       ! neigh_CAD_face = cent_cad_found(bntri2bntri(tet2bn_tri(ii), :))

       ! is_CAD_bn_tri = is_tri_near_CAD_boundary(node2bntri = node2bntri &
       !      , CAD_face = cent_cad_found, nodes = tet_shifted(ii, 1:3))

       do jj = 1, size(rr)
          ! call master2curved_tet( r = r(i), s = s(i), t = t(i), uv = uv &
          !      , xA = xA, x = x(i), y = y(i), z = z(i) )

          ! if ( all(neigh_CAD_face(1) .eq. neigh_CAD_face) ) then

          ! if ( .not.  is_CAD_bn_tri) then

             ! call master2curved_tet_ocas(CAD_face = CAD_face &
             !      , r = rr(jj),s = ss(jj),t = tt(jj) &
             !      , uv = uv, xA = xA, x = xx(jj), y = yy(jj), z = zz(jj))

          ! else

             call master2curved_tet_ocas_close(r = rr(jj),s = ss(jj),t = tt(jj) &
                  , xbot = xbot, xA = xA, tol = tol, x = xx(jj), y = yy(jj), z = zz(jj))

          ! end if

       end do


       ! export the generated curved element
       ! fill this tet coords
       do jj = 1, 4
          xtet(:, jj) = xf(:, tetcon(ii,jj))
       end do
       call tetrahedron_edge_length ( tetra = xtet, edge_length = lens)
       ref_length = 1.5d0 * maxval(lens) / dble(dd)

       if ( indx .eq. 1) then
          call export_tet_face_curve(x = xx, y=yy, z=zz, mina = 0.0d0 &
               , maxa = 160.0d0, fname = 'curved.tec' &
               , meshnum = indx, append_it = .false., ref_length = ref_length)  
       else
          call export_tet_face_curve(x = xx, y=yy, z=zz, mina = 0.0d0 &
               , maxa = 160.0d0, fname = 'curved.tec' &
               , meshnum = indx, append_it = .true., ref_length = ref_length) 
       end if

       indx = indx + 1

       print *, 'indx = ', indx

    end do

    ! map one edge curved tets

    ! map linear tets


    ! clean ups
    if ( allocated(x) ) deallocate(x)
    if ( allocated(icontag) ) deallocate(icontag)
    if ( allocated(xf) ) deallocate(xf)
    if ( allocated(tetcon) ) deallocate(tetcon)
    if ( allocated(neigh) ) deallocate(neigh)
    if ( allocated(bntri) ) deallocate(bntri)

    ! done here
  end subroutine curved_tetgen_geom

  subroutine find_bn_tris_CAD_face(cad_file, nbntri, bntri, xf &
       , cent_cad_found, uvc, tol)
    implicit none
    character(len=*), intent(in) :: cad_file
    integer, intent(in) :: nbntri
    integer, dimension(:), intent(in) :: bntri
    real*8, dimension(:, :), intent(in) :: xf
    integer, allocatable :: cent_cad_found(:) !nbntri
    real*8, allocatable :: uvc(:)
    real*8, intent(in) :: tol

    ! local vars
    ! CAD corresponding data struct
    real*8, allocatable :: xc(:)
    integer :: tpt, ii, jj, kk
    real*8 :: tuv(2), txyz(3), xbn(3,3)

    ! init CAD file
    call init_IGES_f90(fname = cad_file)

    ! find the CAD tag of the centroid of bn faces (tris)
    if ( allocated(cent_cad_found) ) deallocate(cent_cad_found)
    allocate(cent_cad_found(nbntri))
    cent_cad_found = 0

    allocate(xc(3*nbntri))
    xc = 0.0d0

    if ( allocated(uvc) ) deallocate(uvc)
    allocate(uvc(2*nbntri))
    uvc = 0.0d0

    do ii = 1, nbntri
       do jj = 1, 3
          tpt = bntri(6*(ii-1) + jj)
          do kk = 1, 3
             xc(3*(ii-1) + kk) = xc(3*(ii-1) + kk) + xf(kk, tpt)
          end do
       end do
    end do

    ! finalize the center coord
    xc = xc / 3.0d0 
    ! find the CAD tag of the centroids
    print *, 'find the CAD tag of the centroids'
    call find_pts_on_database_f90(npts = nbntri, pts = xc &
         , found = cent_cad_found, uv = uvc, tol = tol)

    ! compute physical center of bn tris and export to MATLAB ...
    open (unit=10, file='tmp.m', status='unknown', action='write')
    write(10, *) 'x = ['
    do ii = 1, nbntri
       tuv(1) = uvc(2*(ii-1) + 1)
       tuv(2) = uvc(2*(ii-1) + 2)
       if (cent_cad_found(ii) .eq. -1) cycle
       call uv2xyz_f90(CAD_face = cent_cad_found(ii), uv = tuv, xyz = txyz)
       print *, 'writing the center of bntri #', ii
       write(10, *) txyz, ';'
    end do
    write(10, *) '];'

    ! print mapped bn triangles
    write(10, *) 'tris = ['
    do ii = 1, nbntri
       if (cent_cad_found(ii) .eq. -1) cycle
       do jj = 1, 3
          tpt = bntri(6*(ii-1) + jj)
          xbn(:, jj) = xf(:, tpt)
       end do
       write(10, *) xbn(:, 1), ';'
       write(10, *) xbn(:, 2), ';'
       write(10, *) xbn(:, 3), ';'
       write(10, *) xbn(:, 1), ';' 
       ! print*, ' '
    end do
    write(10, *) '];'
    close(10)


    ! close CAD objects
    call clean_statics_f90()

    ! cleanups
    if ( allocated(xc) ) deallocate(xc)

    ! done here
  end subroutine find_bn_tris_CAD_face

  !
  subroutine find_tet_adj_csr_zero_based(neigh, xadj, adj)
    implicit none
    integer, dimension(:, :), intent(in) :: neigh
    integer, dimension(:), allocatable :: xadj, adj

    ! local vars
    integer :: ii, jj, nadj, indx

    ! counting and sizing
    nadj = 0
    do ii = 1, size(neigh, 1) ! loop over all tets 
       do jj = 1, size(neigh, 2) ! over neighbors
          if (neigh(ii, jj) > 0 ) nadj = nadj + 1
       end do
    end do

    ! alloc
    if ( allocated(xadj) ) deallocate(xadj)
    allocate(xadj(size(neigh, 1)+1))
    if ( allocated(adj) ) deallocate(adj)
    allocate(adj(nadj))

    ! fill them
    xadj(1) = 0
    indx = 1
    do ii = 1, size(neigh, 1) 
       do jj = 1, size(neigh, 2)
          if (neigh(ii, jj) > 0 ) then
             adj(indx) = neigh(ii, jj)
             indx = indx + 1
          end if
       end do
       xadj(ii+1) = indx - 1
    end do

    adj = adj - 1 ! zero-based cell number

    ! done here
  end subroutine find_tet_adj_csr_zero_based

  !
  subroutine call_metis_graph_parti(xadj, adjncy, vwgt, nparts, part)
    implicit none
    integer, intent(in) :: xadj(:), adjncy(:), vwgt(:)
    integer, intent(in) :: nparts
    integer, intent(out) :: part(:)


    ! local vars
    integer :: nvtxs, ncon
    integer, pointer :: vsize(:) =>null(), adjwgt(:) =>null()
    real*8, pointer :: tpwgts(:)=>null(), ubvec=>null()
    integer, pointer :: options(:) =>null()
    integer :: edgecut

    ! init
    nvtxs = size(xadj) - 1
    ncon = 1

    ! call C-func
    call METIS_PartGraphRecursive(nvtxs, ncon, xadj &
         ,adjncy, vwgt, vsize, adjwgt & 
         ,nparts, tpwgts, ubvec, options & 
         ,edgecut, part)

    ! done here
  end subroutine call_metis_graph_parti

  !
  subroutine shift_tetcon(nbntri, bntri, tetcon, tetcon2, tet2bn_tri)
    implicit none
    integer, intent(in) :: nbntri
    integer, dimension(:), intent(in), target :: bntri
    integer, dimension(:, :), intent(in) :: tetcon
    integer, dimension(:, :), intent(out) :: tetcon2
    integer, dimension(:), intent(out) :: tet2bn_tri

    ! local vars
    integer :: ii, i1, i2, tetnum
    integer, dimension(:), pointer :: pts => null(), tets_on_face => null()

    ! init
    tet2bn_tri = -1

    do ii = 1, nbntri

       i1 = 6* (ii-1) + 1
       i2 = 6* (ii-1) + 3

       pts => bntri(i1:i2)  

       i1 = 6* (ii-1) + 5
       i2 = 6* (ii-1) + 6

       tets_on_face => bntri(i1:i2)

       tetnum = maxval(tets_on_face)

       tet2bn_tri(tetnum) = ii

       ! print *, '***bntri = ', ii, 'bn_pts = ', pts, 'tet_pts = ', tetcon(tetnum, :)
       ! print *, 'results = ', a_in_b(a = pts, b = tetcon(tetnum, 1:3))

       ! bullet proofing ...
       if ( .not. a_in_b(a = pts, b = tetcon(tetnum, :)) ) then
          print *, 'Not all boundary tri points are in the given tet!!! stop'
          stop
       end if

       !
       call shift_tet_to_bn_tri(tet0 = tetcon(tetnum, :), tri = pts &
            , tet = tetcon2(tetnum, :))

!        print *, 'tetcon2 = ', tetcon2(tetnum, :)
! if ( any ( tetcon2(tetnum, :) .ne. tetcon(tetnum, :) ) .and. a_in_b(a = pts, b = tetcon(tetnum, 1:3)) ) then
! print *, 'BADDDDDDDDD'
! stop
! end if

    end do

    ! done here
  end subroutine shift_tetcon

  ! checks see if array(set) <a> is
  ! in array <b>. order is not important.
  ! if "yes" then returns .true. otherwise
  ! returns .false.
  ! 
  function a_in_b(a, b)
    implicit none
    integer, dimension(:), intent(in) :: a, b
    logical :: a_in_b

    ! local vars
    integer :: ii, jj

    do ii = 1, size(a)
       a_in_b = .false.
       do jj = 1, size(b)
          if ( a(ii) .eq. b(jj) ) then
             a_in_b = .true.
             exit
          end if
       end do
       if ( .not. a_in_b ) exit
    end do

    ! done here
  end function a_in_b

  !
  subroutine shift_tet_to_bn_tri(tet0, tri, tet)
    implicit none
    integer, dimension(:), intent(in) :: tet0, tri
    integer, dimension(:), intent(out) :: tet

    ! local vars
    integer :: ii, jj, loc
    logical :: found

    ! init copy
    tet = tet0

    ! first find which tet0(:) point is
    ! not in tri(:). That's the tet's appex.
    ! Then set that as "loc" for shift to right.
    ! The appex should always go to the last
    ! location in connectivity array of the tet
    ! at the end when shifting is complete.

    do ii = 1, 4 

       found = .false.
       do jj = 1, 3
          if ( tri(jj) .eq. tet(ii) ) then
             found = .true.
             exit
          end if
       end do

       if ( .not. found ) then
          loc = ii
          exit
       end if

    end do

    ! Now, shift to right accordingly
    !
    tet = cshift(tet, (ii-4))

    ! done here
  end subroutine shift_tet_to_bn_tri

  subroutine Suv_ocas(uv, CAD_face, S)
    implicit none
    real*8, dimension(:), intent(in) :: uv
    integer, intent(in) :: CAD_face
    real*8, dimension(:), intent(out) :: S

    call uv2xyz_f90(CAD_face = CAD_face, uv = uv, xyz = S)

    ! done here
  end subroutine Suv_ocas

  subroutine master2curved_tet_ocas(CAD_face, r,s,t, uv, xA, x, y, z)
    implicit none
    integer, intent(in) :: CAD_face
    real*8, intent(in) :: r, s, t
    real*8, dimension(:, :), intent(in) :: uv
    real*8, dimension(:), intent(in) :: xA
    real*8, intent(out) :: x, y, z

    ! local vars
    integer :: ii
    real*8 :: u, v, alpha
    real*8, dimension(3) :: Sf, xf
    !
    real*8 :: val   ! the value of basis  
    real*8, dimension(2) :: der, uv_fin 


    if ( abs(t - 1.0d0) <= 1.0d-15 ) then
       u = r ! simple
       v = s ! simple
    else 
       u = r/(1-t)
       v = s/(1-t)
    end if
    alpha = t

    ! compute final uv
    uv_fin = 0.0d0
    do ii = 1, 3
       call psi(etype = 1, i = ii, r = u, s = v, val = val, der = der)
       uv_fin = uv_fin + val * uv(:, ii)
    end do

    ! compute surface points
    ! call Suv(u = uv_fin(1), v= uv_fin(2), S = Sf)
    call Suv_ocas(uv = uv_fin, CAD_face = CAD_face, S = Sf)

    xf = alpha * xA + (1.0d0 - alpha) * Sf

    x = xf(1)
    y = xf(2)
    z = xf(3)

    ! done here
  end subroutine master2curved_tet_ocas

  subroutine master2curved_tet_ocas_close(r,s,t, xbot, xA, tol, x, y, z)
    implicit none
    real*8, intent(in) :: r, s, t
    real*8, dimension(:, :), intent(in) :: xbot
    real*8, dimension(:), intent(in) :: xA
    real*8, intent(in) :: tol
    real*8, intent(out) :: x, y, z

    ! local vars
    integer :: ii, CAD_face(1)
    real*8 :: u, v, alpha
    real*8, dimension(3) :: Sf, xf, xbot_fin
    !
    real*8 :: val   ! the value of basis  
    real*8, dimension(2) :: der, uvout  


    if ( abs(t - 1.0d0) <= 1.0d-15 ) then
       u = r ! simple
       v = s ! simple
    else 
       u = r/(1-t)
       v = s/(1-t)
    end if
    alpha = t

    ! compute final xbot
    xbot_fin = 0.0d0
    do ii = 1, 3
       call psi(etype = 1, i = ii, r = u, s = v, val = val, der = der)
       xbot_fin = xbot_fin + val * xbot(:, ii)
    end do

    ! compute surface points
    ! call Suv(u = uv_fin(1), v= uv_fin(2), S = Sf)
    ! call Suv_ocas(uv = uv_fin, CAD_face = CAD_face, S = Sf)
    call find_pts_on_database_f90(npts = 1, pts = xbot_fin, found = CAD_face, uv = uvout, tol = tol)
    if ( CAD_face(1) .eq. -1 ) then
       print *, 'CAD_face .eq. -1 in master2curved_tet_ocas_close(...)! increase tolerance! stop'
       stop
    end if

    call uv2xyz_f90(CAD_face = CAD_face(1), uv = uvout, xyz = Sf)

    xf = alpha * xA + (1.0d0 - alpha) * Sf

    x = xf(1)
    y = xf(2)
    z = xf(3)

    ! done here
  end subroutine master2curved_tet_ocas_close

  ! subroutine find_bntri2bntri_map(nbntri, bntri, bntri2bntri)
  !   implicit none
  !   integer, intent(in) :: nbntri
  !   integer, dimension(:), intent(in) :: bntri
  !   integer, dimension(:, :), allocatable :: bntri2bntri

  !   ! local vars
  !   integer :: i, i1, i2
  !   integer, allocatable :: bn(:, :) ! bn(1:nbntri, 1:3nodes)

  !   allocate(bn(nbntri, 3))
  !   if (allocated(bntri2bntri) ) deallocate(bntri2bntri)
  !   allocate(bntri2bntri(nbntri, 3))

  !   ! extract "bntri" to a 2d array bn(1:nbntri, 1:3nodes)
  !   ! for ease of work and readability
  !   do i = 1, nbntri
  !      i1 = 6*(i-1) + 1
  !      i2 = 6*(i-1) + 3
  !      bn(i , :) = bntri(i1:i2)
  !   end do

  !   ! now, proceed to fill the output bntri2bntri array
  !   ! using "bn" array info
  !   !
  !   do i = 1, nbntri
  !      bntri2bntri(i, 1) = find_tri_has(bn, i, (/ bn(i, 2), bn(i, 3) /) )
  !      bntri2bntri(i, 2) = find_tri_has(bn, i, (/ bn(i, 3), bn(i, 1) /) )
  !      bntri2bntri(i, 3) = find_tri_has(bn, i, (/ bn(i, 1), bn(i, 2) /) )
  !   end do

  !   ! clean ups
  !   if ( allocated(bn) ) deallocate(bn)

  !   ! done here

  ! contains

  !   function find_tri_has(bn, tri, edg)
  !     implicit none
  !     integer, dimension(:, :) , intent(in) :: bn
  !     integer, intent(in) :: tri
  !     integer, dimension(:), intent(in) :: edg
  !     integer :: find_tri_has

  !     ! local vars
  !     integer :: i, j, cnt

  !     ! init 
  !     find_tri_has = -1 !not found= wall

  !     do i = 1, size(bn, 1)

  !        if ( i .eq. tri ) cycle

  !        cnt = 0
  !        do j = 1, 3
  !           if ( any(bn(i, j) .eq. edg) ) cnt = cnt + 1
  !        end do

  !        if ( cnt .eq. 2) then
  !           find_tri_has = i
  !           exit
  !        end if

  !     end do

  !     ! done here
  !   end function find_tri_has

  ! end subroutine find_bntri2bntri_map

  ! subroutine find_node2bntri_map(nbntri, bntri, node2bntri)
  !   implicit none
  !   integer, intent(in) :: nbntri
  !   integer, dimension(:), intent(in), target :: bntri
  !   type(int_array), dimension(:) :: node2bntri

  !   ! local vars
  !   integer :: i, i1, i2, j
  !   integer, pointer :: nodes(:) => null()

  !   do i = 1, nbntri
  !      i1 = 6*(i-1) + 1
  !      i2 = 6*(i-1) + 3
  !      nodes => bntri(i1:i2)
  !      do j = 1, 3
  !         call push_int_2_array(node2bntri(nodes(j))%val, i)
  !      end do
  !   end do

  !   ! done here

  ! contains

  !   subroutine push_int_2_array(a, i)
  !     implicit none
  !     integer, dimension(:), allocatable :: a
  !     integer, intent(in) :: i

  !     ! local vars
  !     integer :: n
  !     integer, dimension(:), allocatable :: itmp

  !     if ( .not. allocated(a) ) then
  !        n = 1
  !     else
  !        n = size(a) + 1
  !     end if

  !     allocate(itmp(n))
  !     if ( n > 1 ) itmp(1:(n-1)) = a(1:(n-1))
  !     itmp(n) = i
  !     call move_alloc(itmp, a)

  !     ! clean
  !     if ( allocated(itmp) ) deallocate(itmp)

  !     ! done here
  !   end subroutine push_int_2_array

  ! end subroutine find_node2bntri_map

  ! function is_tri_near_CAD_boundary(node2bntri, CAD_face, nodes)
  !   implicit none
  !   type(int_array), dimension(:), intent(in), target :: node2bntri
  !   integer, dimension(:), intent(in) :: CAD_face, nodes
  !   logical :: is_tri_near_CAD_boundary

  !   ! local vars
  !   integer :: inode, ref
  !   integer, pointer :: cells(:) => null()

  !   is_tri_near_CAD_boundary = .false.
  !   ref = CAD_face(node2bntri(nodes(1))%val(1))

  !   do inode = 1, size(nodes)
  !      cells => node2bntri(nodes(inode))%val

  !      if ( .not. all(ref .eq. CAD_face(cells)) ) then
  !         is_tri_near_CAD_boundary = .true.
  !         exit
  !      end if
  !   end do

  !   ! done here
  ! end function is_tri_near_CAD_boundary

end module curved_tet

program tester
  use curved_tet
  implicit none

  ! local vars
  integer :: nhole
  real*8, allocatable :: xh(:)

  !
  ! call tester1()

  ! nhole = 1
  ! allocate(xh(3))
  ! xh = (/ 0.5714d0, 0.4333d0, 0.1180d0 /)

  ! call curved_tetgen_geom(tetgen_cmd = 'pq1.414nnY' &
  !      , facet_file = 'missile_spect3.facet' &
  !      , cad_file = 'store.iges', nhole = nhole, xh = xh, tol = .03d0)

  nhole = 1
  allocate(xh(3))
  xh = (/ 10.0d0, 0.0d0, 0.0d0 /)

  call curved_tetgen_geom(tetgen_cmd = 'pq1.214nnY' &
       , facet_file = 'civil3.facet' &
       , cad_file = 'civil3.iges', nhole = nhole, xh = xh, tol = 20.0d0)

  ! nhole = 1
  ! allocate(xh(3))
  ! xh = 0.0d0

  ! call curved_tetgen_geom(tetgen_cmd = 'pq1.414nnY' &
  !      , facet_file = 'pin.facet' &
  !      , cad_file = 'pin.iges', nhole = nhole, xh = xh, tol = .03d0)

  ! nhole = 1
  ! allocate(xh(3))
  ! xh = 0.0d0

  ! call curved_tetgen_geom(tetgen_cmd = 'pq1.414nnY' &
  !      , facet_file = 'sphere.facet' &
  !      , cad_file = 'sphere2.iges', nhole = nhole, xh = xh, tol = .03d0)

  ! done here
end program tester
