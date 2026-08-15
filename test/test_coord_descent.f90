!
! Correctness test for the coarray coordinate descent.
!
! A synthetic objective with a KNOWN minimum on the simplex, so the answer can
! be checked rather than merely inspected:
!
!   f(w) = sum_k a_k * (w_k - t_k)^2       minimized at w = t
!
! Cheap to evaluate, which is the point -- the parallel line search must be
! validated against a known answer before it is trusted with an objective that
! costs minutes.
!
module quad_obj_mo

  use, intrinsic :: iso_fortran_env, only: real64
  use coord_descent_mo, only: objective_ty

  implicit none

  private
  public :: quad_ty

  integer, parameter :: dp = real64

  type, extends(objective_ty) :: quad_ty
    real(dp), allocatable :: t(:)     ! target (the true optimum)
    real(dp), allocatable :: a(:)     ! per-coordinate curvature
    integer               :: ncall = 0
  contains
    procedure :: evaluate => quad_evaluate
  end type

contains

  subroutine quad_evaluate( this, w, tag, val, ok )
    class(quad_ty), intent(inout) :: this
    real(dp),       intent(in)    :: w(:)
    character(*),   intent(in)    :: tag
    real(dp),       intent(out)   :: val
    logical,        intent(out)   :: ok
    val = sum( this%a * ( w - this%t )**2 )
    this%ncall = this%ncall + 1
    ok = .true.
  end subroutine

end module quad_obj_mo


program test_coord_descent

  use, intrinsic :: iso_fortran_env, only: real64
  use coord_descent_mo, only: cd_ty, cd_result_ty
  use quad_obj_mo,      only: quad_ty

  implicit none

  integer, parameter :: dp = real64
  integer, parameter :: N = 8

  type(cd_ty)        :: cd
  type(quad_ty)      :: obj
  type(cd_result_ty) :: res
  character(64)      :: names(N)
  real(dp) :: w0(N), t(N), a(N), base, err_w, err_f
  integer  :: k
  logical  :: ok
  character(len=8) :: lbl

  do k = 1, N
    write( lbl, '(a,i0)' ) 'coord', k
    names(k) = lbl
    a(k)     = 1.0_dp + 0.3_dp * real( k, dp )    ! distinct curvatures
  end do

  ! A deliberately sparse target: most coordinates want to be at zero, which is
  ! the regime this is used in, and the corner is where a redistribution scheme
  ! is most likely to be wrong.
  t      = 0.0_dp
  t(1)   = 0.50_dp
  t(3)   = 0.30_dp
  t(6)   = 0.20_dp

  ! Start far from the answer: uniform.
  w0 = 1.0_dp / real( N, dp )

  allocate( obj%t, source = t )
  allocate( obj%a, source = a )

  call cd%init( w0, names, delta = 0.05_dp, tol = 1.0e-4_dp, max_pass = 30 )
  cd%verbose = ( this_image() == 1 )

  call obj%evaluate( w0, 'base', base, ok )

  if ( this_image() == 1 ) then
    write( *, '(a)' ) 'coarray coordinate descent -- synthetic quadratic on the simplex'
    write( *, '(a,i0,a,i0)' ) '  coordinates ', N, '   images ', num_images()
    write( *, '(a)' ) '  true optimum:'
    do k = 1, N
      if ( t(k) > 0.0_dp ) write( *, '(a,a10,f10.4)' ) '    ', trim( names(k) ), t(k)
    end do
    write( *, '(a)' ) ''
    flush( 6 )
  end if

  call cd%run( obj, base, res )

  sync all

  if ( this_image() == 1 ) then
    err_w = maxval( abs( res%w - t ) )
    err_f = res%best
    write( *, '(a)' ) ''
    write( *, '(a)' ) '  recovered weights:'
    do k = 1, N
      write( *, '(a,a10,f10.4,a,f10.4)' ) '    ', trim( names(k) ), res%w(k), '   true ', t(k)
    end do
    write( *, '(a,f12.6)' ) '  sum(w)            ', sum( res%w )
    write( *, '(a,f12.6)' ) '  baseline f        ', res%baseline
    write( *, '(a,f12.6)' ) '  final f           ', res%best
    write( *, '(a,f12.6)' ) '  max |w - t|       ', err_w
    write( *, '(a,i0)' )    '  evaluations       ', res%evaluations
    write( *, '(a,i0)' )    '  tolerance-capped  ', res%truncations
    write( *, '(a)' ) ''

    ! The simplex constraint is exact, not approximate -- check it hard.
    if ( abs( sum( res%w ) - 1.0_dp ) > 1.0e-10_dp ) then
      write( *, '(a)' ) '  FAIL  simplex constraint violated'
      error stop 1
    end if
    if ( any( res%w < -1.0e-12_dp ) ) then
      write( *, '(a)' ) '  FAIL  negative weight'
      error stop 1
    end if
    write( *, '(a)' ) '  PASS  sum(w) = 1 exactly, all w >= 0'

    ! delta = 0.05, so the grid cannot resolve better than half a step.
    if ( err_w > 0.05_dp ) then
      write( *, '(a,f8.4,a)' ) '  FAIL  weights off by ', err_w, ' (> one delta)'
      error stop 1
    end if
    write( *, '(a,f8.4,a)' ) '  PASS  weights within one delta of truth (', err_w, ')'

    if ( res%best >= res%baseline ) then
      write( *, '(a)' ) '  FAIL  did not improve on the baseline'
      error stop 1
    end if
    write( *, '(a)' ) '  PASS  improved on the baseline'
    write( *, '(a)' ) ''
    write( *, '(a)' ) '  ALL CHECKS PASSED'
  end if

end program test_coord_descent
