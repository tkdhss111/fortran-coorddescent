!
! User-tunable parameters: each one must actually take effect.
!
! A configuration knob that is silently ignored is worse than one that does not
! exist -- the caller believes they constrained the run and they did not. Each
! case below asserts the OBSERVABLE consequence, not merely that the field was
! accepted.
!
program test_tuning

  use, intrinsic :: iso_fortran_env, only: real64
  use coord_descent_mo, only: cd_ty, cd_result_ty
  use scaled_obj_mo,    only: scaled_ty

  implicit none

  integer, parameter :: dp = real64
  integer, parameter :: N = 8

  integer :: nfail
  nfail = 0

  if ( this_image() == 1 ) write( *, '(a)' ) 'user-tunable parameters'

  call test_order( nfail )
  call test_max_evals( nfail )
  call test_bounds( nfail )
  call test_frozen( nfail )

  if ( this_image() == 1 ) then
    write( *, '(a)' ) ''
    if ( nfail == 0 ) then
      write( *, '(a)' ) '  ALL TUNING CHECKS PASSED'
    else
      write( *, '(a,i0,a)' ) '  ', nfail, ' CHECK(S) FAILED'
      error stop 1
    end if
  end if

contains

  subroutine setup( cd, obj, w0, t )
    type(cd_ty),     intent(inout) :: cd
    type(scaled_ty), intent(inout) :: obj
    real(dp),        intent(out)   :: w0(N), t(N)
    character(64) :: names(N)
    character(8)  :: lbl
    real(dp) :: a(N)
    integer  :: k
    do k = 1, N
      write( lbl, '(a,i0)' ) 'coord', k
      names(k) = lbl
      a(k)     = 1.0_dp + 0.3_dp * real( k, dp )
    end do
    t    = 0.0_dp
    t(1) = 0.50_dp
    t(3) = 0.30_dp
    t(6) = 0.20_dp
    w0   = 1.0_dp / real( N, dp )
    if ( allocated( obj%t ) ) deallocate( obj%t, obj%a )
    allocate( obj%t, source = t )
    allocate( obj%a, source = a )
    obj%scale   = 1.0_dp
    obj%quantum = 0.0_dp
    cd%verbose  = .false.
  end subroutine

  subroutine verdict( label, ok, nf )
    character(*), intent(in)    :: label
    logical,      intent(in)    :: ok
    integer,      intent(inout) :: nf
    if ( this_image() /= 1 ) return
    if ( ok ) then
      write( *, '(a,a)' ) '  PASS  ', label
    else
      write( *, '(a,a)' ) '  FAIL  ', label
      nf = nf + 1
    end if
  end subroutine

  ! Sweep order must reach the search: reversing it must change the trajectory.
  ! Both results stay feasible and converged -- only the path differs.
  subroutine test_order( nf )
    integer, intent(inout) :: nf
    type(cd_ty)     :: cd1, cd2
    type(scaled_ty) :: o1, o2
    type(cd_result_ty) :: r1, r2
    real(dp) :: w0(N), t(N), base
    integer  :: k
    logical  :: ok, differs, feasible

    call setup( cd1, o1, w0, t )
    call cd1%init( w0, [ ( 'c', k = 1, N ) ], delta = 0.05_dp, max_pass = 30 )
    cd1%verbose = .false.
    call o1%evaluate( w0, 'b', base, ok )
    call cd1%run( o1, base, r1 )

    call setup( cd2, o2, w0, t )
    allocate( cd2%order(N) )
    do k = 1, N
      cd2%order(k) = N - k + 1          ! reversed
    end do
    call cd2%init( w0, [ ( 'c', k = 1, N ) ], delta = 0.05_dp, max_pass = 30 )
    cd2%verbose = .false.
    call o2%evaluate( w0, 'b', base, ok )
    call cd2%run( o2, base, r2 )

    sync all
    differs  = any( abs( r1%w - r2%w ) > 1.0e-9_dp )
    feasible = abs( sum( r2%w ) - 1.0_dp ) < 1.0e-9_dp .and. all( r2%w >= -1.0e-12_dp )
    if ( this_image() == 1 ) then
      write( *, '(a,f10.6,a,f10.6)' ) '    forward final ', r1%best, '   reversed final ', r2%best
    end if
    call verdict( 'order changes the trajectory', differs, nf )
    call verdict( 'reversed order stays feasible', feasible, nf )
  end subroutine

  ! A hard budget must be honoured -- that is the whole point of having one.
  subroutine test_max_evals( nf )
    integer, intent(inout) :: nf
    type(cd_ty)     :: cd
    type(scaled_ty) :: obj
    type(cd_result_ty) :: res
    real(dp) :: w0(N), t(N), base
    integer  :: k
    logical  :: ok

    call setup( cd, obj, w0, t )
    cd%max_evals = 30
    call cd%init( w0, [ ( 'c', k = 1, N ) ], delta = 0.05_dp, max_pass = 30 )
    cd%verbose = .false.
    call obj%evaluate( w0, 'b', base, ok )
    call cd%run( obj, base, res )

    sync all
    if ( this_image() == 1 ) write( *, '(a,i0,a)' ) '    evaluations ', res%evaluations, ' (budget 30)'
    ! One block may complete after the budget is reached, so allow a margin of
    ! one full block rather than asserting an exact cut.
    call verdict( 'max_evals bounds the run', res%evaluations <= 30 + 2 * num_images(), nf )
  end subroutine

  ! Per-coordinate bounds must be respected at the final vector.
  subroutine test_bounds( nf )
    integer, intent(inout) :: nf
    type(cd_ty)     :: cd
    type(scaled_ty) :: obj
    type(cd_result_ty) :: res
    real(dp) :: w0(N), t(N), base
    integer  :: k
    logical  :: ok

    call setup( cd, obj, w0, t )
    allocate( cd%wlo(N), source = 0.0_dp )
    allocate( cd%whi(N), source = 1.0_dp )
    cd%wlo(3) = 0.10_dp        ! coord3 may not fall below 0.10
    cd%whi(1) = 0.30_dp        ! coord1 may not exceed 0.30 (true optimum is 0.50)
    call cd%init( w0, [ ( 'c', k = 1, N ) ], delta = 0.05_dp, max_pass = 30 )
    cd%verbose = .false.
    call obj%evaluate( w0, 'b', base, ok )
    call cd%run( obj, base, res )

    sync all
    if ( this_image() == 1 ) then
      write( *, '(a,f8.4,a,f8.4)' ) '    coord1 ', res%w(1), ' (cap 0.30)   coord3 ', res%w(3)
    end if
    call verdict( 'upper bound respected', res%w(1) <= 0.30_dp + 1.0e-9_dp, nf )
  end subroutine

  ! A frozen coordinate must keep its initial weight exactly.
  subroutine test_frozen( nf )
    integer, intent(inout) :: nf
    type(cd_ty)     :: cd
    type(scaled_ty) :: obj
    type(cd_result_ty) :: res
    real(dp) :: w0(N), t(N), base
    integer  :: k
    logical  :: ok

    call setup( cd, obj, w0, t )
    allocate( cd%frozen(N), source = .false. )
    cd%frozen(5) = .true.
    call cd%init( w0, [ ( 'c', k = 1, N ) ], delta = 0.05_dp, max_pass = 30 )
    cd%verbose = .false.
    call obj%evaluate( w0, 'b', base, ok )
    call cd%run( obj, base, res )

    sync all
    if ( this_image() == 1 ) then
      write( *, '(a,f10.6,a,f10.6)' ) '    coord5 initial ', w0(5), '   final ', res%w(5)
    end if
    call verdict( 'frozen coordinate unchanged', abs( res%w(5) - w0(5) ) < 1.0e-9_dp, nf )
    call verdict( 'frozen run stays feasible', abs( sum( res%w ) - 1.0_dp ) < 1.0e-9_dp, nf )
  end subroutine

end program test_tuning
