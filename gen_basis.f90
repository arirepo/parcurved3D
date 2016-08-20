module gen_basis
  ! generates arbitrary 2D basis functions
  ! using Vandermonde matrix operations.
  implicit none

  private

  ! enumerate types for element names grd%elname
  !
  integer, parameter, public :: GEN_TRIANGLE = 1, GEN_QUADRI = 2

  type basis
     private

     real*8, dimension(:,:), allocatable :: MT, M0, U, VT
     integer, dimension(:), allocatable :: piv
     integer :: d
     integer :: elname
     procedure(gen_eval_M_pt), private, pointer, nopass  :: eval_M_pt
     real*8, dimension(:), allocatable :: S
     ! the following holds a copy of original nodes
     ! which is essential for radial basis function operations
     real*8, dimension(:, :), allocatable :: x ! x(1, :) = x, x(2, :) = y
     logical :: is_radial = .false.
     real*8 :: epsil

   contains

     procedure, public :: init => initialize
     procedure, public :: eval => evaluate 
     procedure, private :: svd => comp_svd 
     procedure, public :: fourier_amps => comp_fourier_amps
     procedure, public :: dealloc_basis

  end type basis
  
  abstract interface
     subroutine gen_eval_M_pt(tbasis, xx, yy, d, op, MT)
       import
       implicit none
       class(basis), intent(in) :: tbasis
       real*8, intent(in) :: xx, yy
       integer, intent(in) :: d, op  
       real*8, dimension(:), intent(out) :: MT
     end subroutine gen_eval_M_pt
  end interface

  public :: basis

contains

  ! computes terms in pascal triangle and also their derivatives
  ! with respect to x, y.
  !
  ! inputs:
  !
  ! the terms are evaluated at PHYSICAL point (xx, yy) which can
  ! be inside master element as a special case.
  ! d is the degree of the 2d basis functions  
  ! the operation: (0=pascal), (1=dpascal/dx), (2=dpascal/dy)
  !
  ! outputs:
  !
  ! MT is a real vector containing all pascal terms or
  ! their derivatives d/dx or d/dy at point (xx, yy)
  !
  ! MT = [sum_{i=0}^d sum_{j=0}^i xx**(i-j) * yy**j]
  !
  subroutine comp_pascal_tri(tbasis, xx, yy, d, op, MT)
    implicit none
    class(basis), intent(in) :: tbasis
    real*8, intent(in) :: xx, yy
    integer, intent(in) :: d, op  
    real*8, dimension(:), intent(out) :: MT

    ! local vars
    integer :: id, i, j, jj

    id = size(MT) !

    ! bug checking
    if ( id .ne. ( (d+1)*(d+2)/2 ) ) then
       print *, 'the length of MT does not match with' &
            , ' order of requested polynomial. stop.'
       stop
    end if
    if ( d < 1 ) then
       print *, 'the degree of requested polynomials should be >= 1'
       stop
    end if

    ! fill it!    
    jj = 1
    do i = 0, d
       do j = 0, i

          if     ( op .eq. 0 ) then !pascal terms
             MT(jj) = (xx**(i-j)) * (yy**j)  
          elseif ( op .eq. 1 ) then !d/dx of pascal terms
             if ((i-j-1) < 0) then
                MT(jj) = 0.0d0
             else
                MT(jj) = dble(i-j) * (xx**(i-j-1)) * (yy**j) 
             end if
          elseif ( op .eq. 2 ) then !d/dy of pascal terms
             if ( (j-1) < 0) then
                MT(jj) = 0.0d0
             else
                MT(jj) = (xx**(i-j)) * dble(j) * (yy**(j-1)) 
             end if
          else
             print *, 'unknown operation in comp_pascal_tri! stop.'
             stop 
          end if

          jj = jj + 1
       end do
    end do

    ! done here
  end subroutine comp_pascal_tri

  ! computes basic polynomials for quad elem and also their derivatives
  ! with respect to x, y.
  !
  ! inputs:
  !
  ! the terms are evaluated at PHYSICAL point (xx, yy) which can
  ! be inside master element as a special case.
  ! d is the degree of the 2d basis functions  
  ! the operation: (0=poly), (1=dpoly/dx), (2=dpoly/dy)
  !
  ! outputs:
  !
  ! MT is a real vector containing all polynomial terms or
  ! their derivatives d/dx or d/dy at point (xx, yy)
  !
  ! MT = [sum_{i=0}^d sum_{j=0}^d xx**i * yy**j]
  !
  subroutine comp_poly_quad(tbasis, xx, yy, d, op, MT)
    implicit none
    class(basis), intent(in) :: tbasis
    real*8, intent(in) :: xx, yy
    integer, intent(in) :: d, op  
    real*8, dimension(:), intent(out) :: MT

    ! local vars
    integer :: id, i, j, jj

    id = size(MT) !

    ! bug checking
    if ( id .ne. ( (d+1)**2 ) ) then
       print *, 'the length of MT does not match with' &
            , ' order of requested polynomials for quad! stop.'
       stop
    end if
    if ( d < 1 ) then
       print *, 'the degree of requested polynomials should be >= 1 for quad!'
       stop
    end if

    ! fill it!    
    jj = 1
    do i = 0, d
       do j = 0, d

          if     ( op .eq. 0 ) then !poly terms
             MT(jj) = (xx**i) * (yy**j)  
          elseif ( op .eq. 1 ) then !d/dx of poly
             if (i .eq. 0) then
                MT(jj) = 0.0d0
             else
                MT(jj) = dble(i) * (xx**(i-1)) * (yy**j) 
             end if
          elseif ( op .eq. 2 ) then !d/dy of poly
             if ( j .eq. 0) then
                MT(jj) = 0.0d0
             else
                MT(jj) = (xx**i) * dble(j) * (yy**(j-1)) 
             end if
          else
             print *, 'unknown operation in comp_poly_quad! stop.'
             stop 
          end if

          jj = jj + 1
       end do
    end do

    ! done here
  end subroutine comp_poly_quad

  ! --------------------------------------------
  ! NOTE : The matrix "M" and its transpose "MT" 
  ! used below, is in fact the matrix "V" used 
  ! in module "approx_fekete". It is just different
  ! symbols for the same thing! 
  ! --------------------------------------------
  ! creates MT matrix and then
  ! performs LU factorization of MT matrix.
  ! later this will be used to solve the system and
  ! evaluate basis functions and their 
  ! derivatives at a point 
  subroutine comp_lu_MT(this, x, y)
    implicit none
    class(basis), target, intent(inout) :: this
    real*8, dimension(:), intent(in) :: x, y

    ! local vars
    integer :: i, id, INFO, d, op
    real*8, dimension(:,:), pointer :: MT => null()
    integer, dimension(:), pointer :: piv => null()

    MT  => this%MT
    piv => this%piv
    d = this%d
    op = 0 ! only interpolation

    id = size(x) ! number of points

    ! bug checking
    if ( (id .ne. size(MT,1)) .or. & 
         (id .ne. size(MT,2)) ) then
       print *, 'MT matrix is not (numberofpoints*numberofpoints). stop'
       stop
    end if
    if ( size(piv) .ne. id ) then
       print *, 'the length of pivot vector is not numberofpoints. stop'
       stop
    end if

    ! start filling MT column wise
    do i = 1, id
       call this%eval_M_pt(this, x(i), y(i), d, op, MT(:,i))
    end do

    ! compute SVD and store
    call this%svd()
    ! print *, 'min-max singular value(s) = ', minval(this%S), ' --- ', maxval(this%S)

    ! now, perform lu of MT and store in place
    call DGETRF( id, id, MT, id, piv, INFO )
    if (INFO .ne. 0) then
       print *, 'somethig is wrong in LU decomposition in MT! stop.'
       stop
    end if

    ! done here
  end subroutine comp_lu_MT

  
  ! initializes the basis data type
  subroutine initialize(this, x, y, elname, do_radial)
    implicit none
    class(basis), intent(inout) :: this
    real*8, dimension(:), intent(in) :: x, y
    integer, intent(in) :: elname
    logical, intent(in), optional :: do_radial

    ! local vars
    integer :: id
    real*8  :: delta

    id = size(x)

    ! bulletproofing
    if ( allocated(this%MT) ) then
       print *, 'this object for basis function is already initialized! stop'
       stop
    end if

    !
    allocate(this%MT(id, id), this%M0(id, 1), this%piv(id))
    this%M0 = 0.0d0
    this%elname = elname
    this%is_radial = .false.

    ! take a copy of the input nodes 
    allocate(this%x(2, id))
    this%x(1, :) = x
    this%x(2, :) = y
        
    ! obtaining the degree of required polynomial
    ! from the given number of points and performing
    ! element specific initialization
    select case (this%elname)
    case (GEN_TRIANGLE ) ! do it the old way

       delta = sqrt(1.0d0 + 8.0d0 * dble(id))
       this%d = maxval( (/ nint(-1.5d0 + 0.5d0 * delta) &
            , nint(-1.5d0 - 0.5d0 * delta) /) )

       if ( this%d <= 0 ) then
          print *, 'degree of basis is this%d <= 0! stop.'
          stop
       end if

       ! initialize generic evaluation procedure pointer
       this%eval_M_pt => comp_pascal_tri

    case ( GEN_QUADRI )

       this%d = nint(sqrt(dble(id))) - 1

       ! initialize generic evaluation procedure pointer
       this%eval_M_pt => comp_poly_quad

    case default

       print *, 'unregognized element name in init. basis func.! stop'
       stop

    end select

    ! overwrite with radial basis settings if radial option is selected
    if ( present( do_radial ) ) then
       if ( do_radial ) then !really do radial!

          this%d = id
          ! initialize generic evaluation procedure pointer
          this%eval_M_pt => comp_radial_terms
          ! set the radial interpolant property
          this%is_radial = .true.
          ! find approximately best RBF parameter
          call comp_RBF_epsil(this)

       end if
    end if

    ! fill MT matrice
    call comp_lu_MT(this, x, y)

    ! done here
  end subroutine initialize

  ! evaluates the basis function or
  ! it derivatives d/dx, d/dy at 
  ! some point (x0, y0)
  ! op = 0 => evaluate the basis function
  ! op = 1 => evaluate the d/dx of basis
  ! op = 2 => evaluate the d/dy of basis
  subroutine evaluate(this, x0, y0, op, val, resol)
    implicit none
    class(basis), intent(inout) :: this
    real*8, intent(in) :: x0, y0
    integer, intent(in) :: op
    real*8, dimension(:), intent(out) :: val
    real*8, intent(in), optional :: resol

    ! local vars
    integer :: id, INFO, km

    ! bullet proofing
    if (.not. allocated(this%MT)) then
       print *, 'fatal: please first initialize this basis object' &
            , ' before evaluating. stop.'
       stop
    end if

    id = size(this%MT, 1)

    if ( size(val) .ne. id) then
       print *, 'in evaluating basis functions, size of' &
            , ' output array size<val> = ', size(val), ' is not equal to' &
            , ' the number of basis functions (numberofpoints = ', id, '). stop.'
       stop
    end if

    ! first fill this%M0 column vector with
    ! pascal terms at point (x0, y0) according to the requested operation
    call this%eval_M_pt(this, x0, y0, this%d, op, val)

    if ( .not. present(resol) ) then
       ! solve psi = basis = MT\M0 using already computed LU
       CALL DGETRS( 'No transpose', id, 1, this%MT, id, this%piv, val, id, INFO )
       if ( INFO .ne. 0) then
          print *, 'could not solve to evaluate basis at point (x0, y0)'
          stop
       end if
    else

       km = nint( dble(id) * resol)
       if ( (km <= 0) .or. (km > id) ) then
          print *, 'fatal : wrong maximum Fourier mode number! stop'
          stop
       end if

       val = matmul(this%U(:, 1:km) , (matmul(this%VT(1:km, :), val) / this%S(1:km)))

    end if

    ! val = this%M0(:,1)

    ! done here
  end subroutine evaluate

  ! computes the SVD of M = transpose(MT) matrix 
  ! which will be used for modal analysis
  ! and spectral filtering
  subroutine comp_svd(this)
    implicit none
    class(basis), intent(inout) :: this

    ! local vars
    character :: JOBZ
    integer :: M, N, LDA, LDU, LDVT, LWORK, INFO
    real*8, dimension(:, :), allocatable :: A
    real*8, dimension(:), allocatable :: WORK
    integer, dimension(:), allocatable :: IWORK

    ! bullet proofing
    if (.not. allocated(this%MT)) then
       print *, 'MT needs to be init. before SVD! stop.'
       stop
    end if

    ! LAPACK signature for gen. SVD
    !
    ! subroutine dgesdd(character JOBZ,
    !   integer M,
    !   integer N,
    !   double precision, dimension( lda, * ) A,
    !   integer LDA,
    !   double precision, dimension( * ) S,
    !   double precision, dimension( ldu, * ) U,
    !   integer LDU,
    !   double precision, dimension( ldvt, * ) VT,
    !   integer LDVT,
    !   double precision, dimension( * ) WORK,
    !   integer LWORK,
    !   integer, dimension( * ) IWORK,
    !   integer INFO 
    !   )

    ! init.
    allocate( A( size(this%MT,2), size(this%MT,1) ) )
    A = transpose(this%MT)

    JOBZ = 'A'
    M = size(A, 1)
    N = size(A, 2)
    LDA = max(1, M)
    LDU = M
    LDVT = N
    LWORK = 2 * ( min(M,N)*(6+4*min(M,N))+max(M,N) )


    if ( allocated(this%S) ) deallocate(this%S)
    allocate(this%S(min(M, N)))

    if ( allocated(this%U) ) deallocate(this%U)
    allocate(this%U(LDU, M))

    if ( allocated(this%VT) ) deallocate(this%VT)
    allocate(this%VT(LDVT, N))

    allocate(WORK(max(1, LWORK)))
    allocate(IWORK( (8*min(M, N)) ))


    ! compute and store SVD in this object
    call dgesdd(JOBZ, M, N, A, LDA, this%S &
         ,this%U, LDU, this%VT, LDVT, WORK, LWORK &
         , IWORK, INFO)

    if (INFO .ne. 0) then
       print *, 'something is wrong in SVD decomposition of MT! INFO = ', INFO, ' . stop'
       stop
    end if

    ! clean ups
    deallocate(A, WORK, IWORK)

    ! done here
  end subroutine comp_svd

  subroutine comp_fourier_amps(this, u0, amps)
    implicit none
    class(basis), intent(in) :: this
    real*8, dimension(:), intent(in) :: u0
    real*8, dimension(:), intent(out) :: amps

    amps = matmul(transpose(this%U), u0)

    ! done here
  end subroutine comp_fourier_amps
  
  !
  ! generates the radial terms for "d" number
  ! of scattered points.
  !
  ! NOTE : here "d" is not the order of approximation
  ! and is just the number of points.
  !
  subroutine comp_radial_terms(tbasis, xx, yy, d, op, MT)
    implicit none
    class(basis), intent(in) :: tbasis
    real*8, intent(in) :: xx, yy
    integer, intent(in) :: d, op  
    real*8, dimension(:), intent(out) :: MT

    ! local vars
    integer :: i
    real*8 :: r, xi, yi, u

    ! bullet proofing
    if ( d .ne. size(tbasis%x , 2) ) then
       print *, ' the given parameter <d> is not equal to the' &
            , ' number of scattered points! stop'
       stop
    end if

    ! compute and fill terms 
    do i = 1, d

       xi = tbasis%x(1, i)
       yi = tbasis%x(2, i)

       r = sqrt( ( xx - xi )**2 + ( yy - yi )**2 )
       u = exp( -1.0d0 * (tbasis%epsil * r)**2 )

       if ( op .eq. 0 ) then !interp.
          MT(i) = u
       elseif (op .eq. 1) then !d/dx
          MT(i) = -2.0d0 * (tbasis%epsil**2) * (xx - xi) * u
       elseif (op .eq. 2) then !d/dy
          MT(i) = -2.0d0 * (tbasis%epsil**2) * (yy - yi) * u
       else
          print *, 'unknown operation in comp_radial_terms! stop.'
          stop 
       end if

    end do

    ! done here
  end subroutine comp_radial_terms

  ! finds the root of function f(x) in interval [a, b]
  ! satisfying |f(z)| <= eps_abs.
  ! The result is returned in "r".
  subroutine bisection(this, f, a, b, eps_abs, r )
    implicit none
    class(basis) :: this
    real*8 :: f, a, b
    real*8, intent(in) :: eps_abs
    real*8, intent(out) :: r

    ! local vars
    real*8 :: c, fa, fb, fc

    ! bullet proofing ...
    ! Check that that neither end-point is a root
    ! and if f(a) and f(b) have the same sign, throw an exception.
    if ( abs(f(this,a)) <= eps_abs ) then
       r = a
       return
    elseif ( abs(f(this,b)) <= eps_abs ) then
       r = b
       return
    elseif ( (f(this,a) * f(this,b)) > 0.0d0 ) then
       print *, 'f(a) and f(b) do not have opposite signs! stop'
       stop
    end if

    ! We will iterate until bisection converges!
    ! NOTE : always choose rough tolerances to increase
    ! the efficiency 
    do 

       ! Find the mid-point
       c = (a + b)/2.0d0

       ! evaluate ONCE
       fa = f(this,a)
       fb = f(this,b)
       fc = f(this,c)

       ! Check if we found a root or whether or not
       ! we should continue with:
       !          [a, c] if f(a) and f(c) have opposite signs, or
       !          [c, b] if f(c) and f(b) have opposite signs.
       if ( abs(fc) <= eps_abs ) then
          r = c
          return
       elseif ( fc*fa < 0.0d0 ) then
          b = c
       else
          a = c
       end if

       ! check the possibly converged root
       !
       ! |f(a)| < |f(b)| and |f(a)| < eps_abs and return 'a', or
       ! |f(b)| < eps_abs and return 'b'.
       if ( (abs( fa ) < abs( fb )) .and. (abs( fa ) <= eps_abs) ) then
          r = a
          return
       elseif ( abs( fb ) <= eps_abs ) then
          r = b
          return
       end if

    end do !main loop

    ! the following line will not be reached!
    print *, 'Warning: bisection did not converge!'

    ! done here
  end subroutine bisection

  ! dynamically computes the minimum singular value
  ! of the Vandermonde matrix of RBF basis for the
  ! given RBF parameter "epsil"
  function min_sing(this, epsil)
    implicit none
    class(basis), target :: this
    real*8, intent(in) :: epsil
    real*8 :: min_sing

     ! local vars
    integer :: i, id, d, op
    real*8, dimension(:,:), pointer :: MT => null()

    ! init
    MT  => this%MT
    d = this%d
    op = 0 ! only interpolation
    this%epsil = epsil

    id = size(this%x, 2) ! number of points

    ! bug checking
    if ( (id .ne. size(MT,1)) .or. & 
         (id .ne. size(MT,2)) ) then
       print *, 'MT matrix is not (numberofpoints*numberofpoints). stop'
       stop
    end if

    ! start filling MT column wise
    do i = 1, id
       call this%eval_M_pt(this, this%x(1,i), this%x(2,i), d, op, MT(:,i))
    end do

    ! compute SVD and store
    call this%svd()

    ! now return the minimum singular value
    min_sing = minval(this%S)
 
    ! done here
  end function min_sing

  !
  ! TARGET FUNCTION
  !
  ! for finding zeros of equation:
  !   sigma_min = (close to machine eps)
  ! which is used as the approx. best RBF param.
  !
  ! NOTE : we permanently chose 
  ! (close to machine eps) = 1.0d-15
  ! for double precision
  !
  function targ(this, epsil)
    implicit none
    class(basis) :: this
    real*8, intent(in) :: epsil
    real*8 :: targ

    targ = min_sing(this, epsil) - 1.0d-14

    ! done here
  end function targ

  ! computes approximately best RBF param.
  !  
  subroutine comp_RBF_epsil(this)
    implicit none
    class(basis) :: this

    ! local vars
    real*8 :: a, b, eps_abs, r

    a = 0.0001d0 ! small value
    b = dble(size(this%x, 2)) ! nmax
    eps_abs = 1.d-15

    ! find approx. best RBF param, i.e. epsil
    call bisection(this, targ, a, b, eps_abs, r )

    ! save it into this object
    this%epsil = r

    ! done here
  end subroutine comp_RBF_epsil

  ! completely cleans/deallocates a basis
  ! function struct
  !
  subroutine dealloc_basis(tbasis)
    implicit none
    class(basis), intent(inout) :: tbasis


    if(allocated(tbasis%MT)) deallocate(tbasis%MT)
    if(allocated(tbasis%M0)) deallocate(tbasis%M0)
    if(allocated(tbasis%piv)) deallocate(tbasis%piv)

    if(associated(tbasis%eval_M_pt)) nullify(tbasis%eval_M_pt)

    ! done here
  end subroutine dealloc_basis

end module gen_basis

! program tester
!   use gen_basis
!   implicit none

!   integer :: i, npelem
!   real*8 :: x0, y0
!   real*8, dimension(:), allocatable :: xi, eta

!   ! first order element
!   print *, 'first order element :'
!   npelem = 3
!   allocate( xi(npelem), eta(npelem))

!   xi  = (/ 0.0d0, 1.0d0, 0.0d0 /)
!   eta = (/ 0.0d0, 0.0d0, 1.0d0 /)  
  
!   do i = 1, size(xi)

!      x0 = xi(i); y0 = eta(i)
!      call print_basis(xi, eta, x0, y0)

!   end do

!   deallocate(xi, eta)

!   ! second order element
!   print *, 'second order element:'
!   npelem = 6
!   allocate( xi(npelem), eta(npelem))

!   xi  = (/ 0.0d0, 0.5d0, 1.0d0, 0.5d0, 0.0d0, 0.0d0 /)
!   eta = (/ 0.0d0, 0.0d0, 0.0d0, 0.5d0, 1.0d0, 0.5d0 /)  
  
!   do i = 1, size(xi)

!      x0 = xi(i); y0 = eta(i)
!      call print_basis(xi, eta, x0, y0)

!   end do

!   deallocate(xi, eta)

!   ! done here

! contains

!   subroutine print_basis(xi, eta, x0, y0)
!     implicit none
!     real*8, dimension(:), intent(in) :: xi, eta
!     real*8, intent(in) :: x0, y0

!     ! local vars
!     type(basis) :: tbasis
!     real*8, dimension(size(xi)) :: val

!     call tbasis%init(xi, eta)

!     print *, 'at point (x0, y0) = ', x0, ',', y0
!     call tbasis%eval(x0, y0, 0, val)
!     print *, 'psi = ', val
!     call tbasis%eval(x0, y0, 1, val)
!     print *, 'd psi/dx = ', val
!     call tbasis%eval(x0, y0, 2, val)
!     print *, 'd psi/dy = ', val

!     ! done here
!   end subroutine print_basis

! end program tester
