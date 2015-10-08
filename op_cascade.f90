module op_cascade
  implicit none

  private

  interface
     subroutine init_STEP_intf(fname) &
          bind( C, name = "init_STEP")
       use iso_c_binding, only : c_char
       import
       implicit none
       character (c_char), intent(in) :: fname(*)

     end subroutine init_STEP_intf

  end interface

  interface
     subroutine init_IGES_intf(fname) &
          bind( C, name = "init_IGES")
       use iso_c_binding, only : c_char
       import
       implicit none
       character (c_char), intent(in) :: fname(*)

     end subroutine init_IGES_intf

  end interface

  interface
     subroutine clean_statics_intf() &
          bind( C, name = "clean_statics")
       implicit none


     end subroutine clean_statics_intf

  end interface

  interface
     subroutine find_pts_on_database_intf(npts, pts, found, uv, tol) &
          bind( C, name = "find_pts_on_database")
       use iso_c_binding, only : c_int, c_double
       import
       implicit none
       integer(c_int), value :: npts
       real(c_double) :: pts(3 * npts)
       integer(c_int) :: found(npts)
       real(c_double) :: uv(2*npts)
       real(c_double), value :: tol

     end subroutine find_pts_on_database_intf

  end interface

  interface
     subroutine uv2xyz_intf(CAD_face, uv, xyz) &
          bind( C, name = "uv2xyz")
       use iso_c_binding, only : c_int, c_double
       import
       implicit none
       integer(c_int), value :: CAD_face
       real(c_double) :: uv(2)
       real(c_double) :: xyz(3)

     end subroutine uv2xyz_intf

  end interface

  interface
     subroutine xyz2uv_intf(CAD_face, xyz, uv, tol) &
          bind( C, name = "xyz2uv")
       use iso_c_binding, only : c_int, c_double
       import
       implicit none
       integer(c_int), value :: CAD_face
       real(c_double) :: xyz(3)
       real(c_double) :: uv(2)
       real(c_double), value :: tol

     end subroutine xyz2uv_intf

  end interface

  ! Readers
  public :: init_STEP_f90, init_IGES_f90
  ! Queries
  public :: find_pts_on_database_f90, uv2xyz_f90, xyz2uv_f90
  ! Cleanups
  public :: clean_statics_f90

contains

  subroutine init_STEP_f90(fname)
    use iso_c_binding, only : c_char, c_null_char
    implicit none
    character(len = *), intent(in) :: fname

    !local vars
    character (kind = c_char, len = 255) :: fname_in
    fname_in = trim(fname)//C_NULL_CHAR 

    ! call C-func
    call init_STEP_intf(fname_in)

    ! done here
  end subroutine init_STEP_f90

  subroutine init_IGES_f90(fname)
    use iso_c_binding, only : c_char, c_null_char
    implicit none
    character(len = *), intent(in) :: fname

    !local vars
    character (kind = c_char, len = 255) :: fname_in
    fname_in = trim(fname)//C_NULL_CHAR 

    ! call C-func
    call init_IGES_intf(fname_in)

    ! done here
  end subroutine init_IGES_f90

  subroutine clean_statics_f90()
    implicit none

    ! call C-func
    call clean_statics_intf()

    ! done here
  end subroutine clean_statics_f90

  subroutine find_pts_on_database_f90(npts, pts, found, uv, tol)
    use iso_c_binding, only : c_int, c_double
    implicit none
    integer, intent(in) :: npts
    real*8, dimension(:), intent(in) :: pts
    integer, dimension(:), intent(out) :: found
    real*8, dimension(:), intent(out) :: uv
    real*8, intent(in) :: tol

    !local vars
    integer :: i
    integer(c_int) :: npts_in
    real(c_double) :: tol_in
    real(c_double), allocatable :: pts_in(:)
    integer(c_int), allocatable :: found_out(:)
    real(c_double), allocatable :: uv_out(:)

    ! init 
    npts_in = npts
    tol_in = tol

    allocate(pts_in(size(pts)))
    do i = 1, size(pts)
       pts_in(i) = pts(i)
    end do

    allocate(found_out(size(found)))
    found_out = 0
    allocate(uv_out(size(uv)))
    uv_out = 0.0d0

    ! call C-func
    call find_pts_on_database_intf(npts_in, pts_in, found_out, uv_out, tol_in)

    ! fillout the return arrays and values
    do i = 1, size(found)
       found(i) = found_out(i)
    end do

    do i = 1, size(uv)
       uv(i) = uv_out(i)
    end do

    ! clean ups
    if ( allocated(pts_in) ) deallocate(pts_in)
    if ( allocated(found_out) ) deallocate(found_out)
    if ( allocated(uv_out) ) deallocate(uv_out)

    ! done here
  end subroutine find_pts_on_database_f90

  subroutine uv2xyz_f90(CAD_face, uv, xyz)
    use iso_c_binding, only : c_int, c_double
    implicit none
    integer, intent(in) :: CAD_face
    real*8, dimension(:), intent(in) :: uv
    real*8, dimension(:), intent(out) :: xyz

    !local vars
    integer :: i
    integer(c_int) :: CAD_face_in
    real(c_double) :: uv_in(2)
    real(c_double) :: xyz_out(3)

    ! init 
    CAD_face_in = CAD_face

    do i = 1, 2
       uv_in(i) = uv(i)
    end do

    ! call C-func
    call uv2xyz_intf(CAD_face_in, uv_in, xyz_out)

    ! fillout the return arrays and values
    do i = 1, 3
       xyz(i) = xyz_out(i)
    end do

    ! done here
  end subroutine uv2xyz_f90

  subroutine xyz2uv_f90(CAD_face, xyz, uv, tol)
    use iso_c_binding, only : c_int, c_double
    implicit none
    integer, intent(in) :: CAD_face
    real*8, dimension(:), intent(in) :: xyz
    real*8, dimension(:), intent(out) :: uv
    real*8, intent(in) :: tol

    !local vars
    integer :: i
    integer(c_int) :: CAD_face_in
    real(c_double) :: xyz_in(3)
    real(c_double) :: uv_out(2)
    real(c_double) :: tol_in

    ! init 
    CAD_face_in = CAD_face

    do i = 1, 3
       xyz_in(i) = xyz(i)
    end do

    tol_in = tol

    ! call C-func
    call xyz2uv_intf(CAD_face_in, xyz_in, uv_out, tol_in)

    ! fillout the return arrays and values
    do i = 1, 2
       uv(i) = uv_out(i)
    end do

    ! done here
  end subroutine xyz2uv_f90

end module op_cascade

! program tester
!   use op_cascade
!   implicit none

!   ! local vars
!   character(len = 255) :: fname_step = './crude_samples/store.step'
!   character(len = 255) :: fname_iges = './crude_samples/store.iges'
!   ! chose a point in 3D space for query
!   integer :: npts
!   real*8 :: pts(3)
!   integer :: found(1)
!   real*8 :: uv(2)
!   real*8 :: tol

!   npts = 1
!   pts = (/0.373333d0, 0.433333d0, 0.118d0 /)
!   found = 0
!   uv = 0.0d0
!   tol = 1.0d-1

!   ! read STEP file and init.
!   call init_STEP_f90(fname_step)

!   ! do a query
!   call find_pts_on_database_f90(npts, pts, found, uv, tol)

!   ! show result
!   print *, ' ===========  in STEP file: =========== '
!   print *, 'found(1) = ', found(1)
!   print *, '(u, v) = (', uv(1), ', ',  uv(2), ') '

!   ! Now, test snapping function uv2xyz
!   pts = (/ 0.0d0, 0.0d0, 0.0d0 /)
!   call uv2xyz_f90(found(1), uv, pts)
!   print *, 'snapped point = ', pts
 
!   ! final cleanups
!   call clean_statics_f90()

!   ! done here
! end program tester
