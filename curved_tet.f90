module curved_tet
  use tetmesher
  use tet_props
  use lag_basis
  use op_cascade
  implicit none

  private

  real*8, parameter :: Radius = 2.0d0

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
  subroutine export_tet_face_curve(x, y, z, mina, maxa, fname, meshnum, append_it)
    implicit none
    real*8, dimension(:), intent(in) :: x, y, z
    real*8, intent(in) :: mina, maxa
    character(len = *), intent(in) :: fname
    integer, intent(in) :: meshnum
    logical , intent(in) :: append_it

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
    call filter_bad_tets(x = xf, icon = tetcon, mina = mina &
         , maxa= maxa, active = is_active)

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
  subroutine filter_bad_tets(x, icon, mina, maxa, active)
    implicit none
    real*8, dimension(:, :), intent(in) :: x
    integer, dimension(:, :), intent(in) :: icon
    real*8, intent(in) :: mina, maxa
    logical, dimension(:), intent(out) :: active

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
    integer :: npts, nquad, ntri
    real*8, dimension(:), allocatable :: x
    integer, dimension(:), allocatable :: icontag

    ! outs
    real*8, dimension(:, :), allocatable :: xf
    integer, dimension(:, :), allocatable :: tetcon, neigh
    integer :: nbntri
    integer, dimension(:), allocatable :: bntri

    real*8, allocatable :: uu(:, :)

    ! CAD corresponding data struct
    integer, allocatable :: cent_cad_found(:) !nbntri
    real*8, allocatable :: xc(:), uvc(:)
    integer :: tpt, ii, jj, kk
    real*8 :: tuv(2), txyz(3), xbn(3,3)

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

    ! export linear tetmesh
    allocate(uu(1, size(xf,2)))
    uu = 1.0d0
    call write_u_tecplot_tet(meshnum=1, outfile='linear_tets.tec', x = xf &
         , icon = tetcon, u = uu, appendit = .false.)
    if ( allocated(uu) ) deallocate(uu)

    ! init CAD file
    call init_IGES_f90(fname = cad_file)

    ! find the CAD tag of the centroid of bn faces (tris)
    allocate(cent_cad_found(nbntri))
    cent_cad_found = 0
    allocate(xc(3*nbntri))
    xc = 0.0d0
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

    ! print *, 'cent_cad_found = ', cent_cad_found
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

    ! map one face curved tets

    ! map one edge curved tets

    ! close CAD objects
    call clean_statics_f90()

    ! clean ups
    if ( allocated(x) ) deallocate(x)
    if ( allocated(icontag) ) deallocate(icontag)
    if ( allocated(xf) ) deallocate(xf)
    if ( allocated(tetcon) ) deallocate(tetcon)
    if ( allocated(neigh) ) deallocate(neigh)
    if ( allocated(bntri) ) deallocate(bntri)
    if ( allocated(cent_cad_found) ) deallocate(cent_cad_found)

    ! done here
  end subroutine curved_tetgen_geom

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

  call curved_tetgen_geom(tetgen_cmd = 'pq1.414nnY' &
       , facet_file = 'civil3.facet' &
       , cad_file = 'civil3.iges', nhole = nhole, xh = xh, tol = 20.0d0)

  ! done here
end program tester
