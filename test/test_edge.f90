!
! Edge cases.
!
! These are chosen for the failure modes that do not produce a wrong number --
! they hang, crash, or loop forever. A search that returns a bad answer is
! recoverable; one that never returns is not, and one that crashes under
! coarrays returns exit code 0 while doing it.
!
! Every case must TERMINATE and leave a feasible vector. That is the bar.
!
module edge_obj_mo

  use, intrinsic :: iso_fortran_env, only: real64
  use coord_descent_mo, only: objective_ty

  implicit none

  private
  public :: edge_ty

  integer, parameter :: dp = real64

  ! kind: 1 quadratic, 2 constant, 3 NaN, 4 always-fail, 5 optimum outside simplex
  type, extends(objective_ty) :: edge_ty
    integer               :: kind  = 1
    integer               :: ncall = 0
    real(dp), allocatable :: t(:)
  contains
    procedure :: evaluate => edge_evaluate
  end type

contains

  subroutine edge_evaluate( this, w, tag, val, ok )
    class(edge_ty), intent(inout) :: this
    real(dp),       intent(in)    :: w(:)
    character(*),   intent(in)    :: tag
    real(dp),       intent(out)   :: val
    logical,        intent(out)   :: ok
    real(dp) :: zero

    this%ncall = this%ncall + 1
    ok  = .true.
    val = 0.0_dp

    select case ( this%kind )
    case ( 1 )
      val = sum( ( w - this%t )**2 )
    case ( 2 )
      val = 1.0_dp                       ! constant: no gradient anywhere
    case ( 3 )
      zero = 0.0_dp
      val  = zero / zero                 ! NaN
    case ( 4 )
      val = 0.0_dp
      ok  = .false.                      ! every evaluation fails
    case ( 5 )
      val = sum( ( w - 5.0_dp )**2 )     ! optimum far outside the simplex
    end select
  end subroutine

end module edge_obj_mo


program test_edge

  use, intrinsic :: iso_fortran_env, only: real64
  use coord_descent_mo, only: cd_ty, cd_result_ty
  use edge_obj_mo,      only: edge_ty

  implicit none

  integer, parameter :: dp = real64
  integer :: nfail
  nfail = 0

  if ( this_image() == 1 ) write( *, '(a)' ) 'edge cases'

  call case_min_dim( nfail )
  call case_start_at_optimum( nfail )
  call case_start_at_corner( nfail )
  call case_constant_objective( nfail )
  call case_nan_objective( nfail )
  call case_all_fail( nfail )
  call case_huge_delta( nfail )
  call case_tiny_delta( nfail )
  call case_one_pass( nfail )
  call case_unnormalized_input( nfail )
  call case_optimum_outside( nfail )

  if ( this_image() == 1 ) then
    write( *, '(a)' ) ''
    if ( nfail == 0 ) then
      write( *, '(a)' ) '  ALL EDGE CASES PASSED'
    else
      write( *, '(a,i0,a)' ) '  ', nfail, ' EDGE CASE(S) FAILED'
      error stop 1
    end if
  end if

contains

  subroutine verdict( label, ok, nf, extra )
    character(*), intent(in)    :: label
    logical,      intent(in)    :: ok
    integer,      intent(inout) :: nf
    character(*), optional, intent(in) :: extra
    if ( this_image() /= 1 ) return
    if ( ok ) then
      write( *, '(a,a40,a)' ) '  PASS  ', label, merge_txt( extra )
    else
      write( *, '(a,a40,a)' ) '  FAIL  ', label, merge_txt( extra )
      nf = nf + 1
    end if
  end subroutine

  function merge_txt( extra ) result( t )
    character(*), optional, intent(in) :: extra
    character(80) :: t
    t = ''
    if ( present( extra ) ) t = extra
  end function

  logical function feasible( w )
    real(dp), intent(in) :: w(:)
    feasible = abs( sum( w ) - 1.0_dp ) < 1.0e-9_dp .and. all( w >= -1.0e-12_dp ) &
               .and. all( w == w )
  end function

  ! Runs one configuration and returns the result; every case funnels through
  ! here so "did it terminate" is uniform.
  subroutine run_case( n, kind, w0, delta, max_pass, res, ncall )
    integer,  intent(in)  :: n, kind, max_pass
    real(dp), intent(in)  :: w0(:), delta
    type(cd_result_ty), intent(out) :: res
    integer,  intent(out) :: ncall

    type(cd_ty)   :: cd
    type(edge_ty) :: obj
    character(64), allocatable :: names(:)
    real(dp) :: base
    integer  :: k
    logical  :: ok

    allocate( names(n) )
    do k = 1, n
      write( names(k), '(a,i0)' ) 'c', k
    end do
    obj%kind = kind
    ! The target must lie ON the simplex, or "start at the optimum" starts
    ! somewhere else and the case tests nothing it claims to.
    allocate( obj%t(n) )
    obj%t(1)   = 0.4_dp
    obj%t(2:n) = 0.6_dp / real( n - 1, dp )

    call cd%init( w0, names, delta = delta, max_pass = max_pass )
    cd%verbose = .false.
    call obj%evaluate( w0, 'b', base, ok )
    call cd%run( obj, base, res )
    sync all
    ncall = obj%ncall
  end subroutine

  ! Two free coordinates is the minimum the search can act on.
  subroutine case_min_dim( nf )
    integer, intent(inout) :: nf
    type(cd_result_ty) :: res
    integer :: nc
    call run_case( 2, 1, [ 0.5_dp, 0.5_dp ], 0.05_dp, 10, res, nc )
    call verdict( 'n = 2 (minimum dimension)', feasible( res%w ), nf )
  end subroutine

  ! Starting AT the optimum must converge immediately, not wander.
  subroutine case_start_at_optimum( nf )
    integer, intent(inout) :: nf
    type(cd_result_ty) :: res
    real(dp) :: w0(4)
    integer :: nc
    w0      = 0.6_dp / 3.0_dp
    w0(1)   = 0.4_dp                     ! exactly the target used in run_case
    call run_case( 4, 1, w0, 0.05_dp, 10, res, nc )
    call verdict( 'start at the optimum', feasible( res%w ) .and. res%best <= 1.0e-6_dp, nf )
  end subroutine

  ! A corner of the simplex: one coordinate holds everything, the rest are 0.
  ! This is where redistribution has an empty pool to divide.
  subroutine case_start_at_corner( nf )
    integer, intent(inout) :: nf
    type(cd_result_ty) :: res
    real(dp) :: w0(5)
    integer :: nc
    w0    = 0.0_dp
    w0(1) = 1.0_dp
    call run_case( 5, 1, w0, 0.05_dp, 10, res, nc )
    call verdict( 'start at a simplex corner', feasible( res%w ), nf )
  end subroutine

  ! No gradient anywhere. Must terminate rather than sweep forever.
  subroutine case_constant_objective( nf )
    integer, intent(inout) :: nf
    type(cd_result_ty) :: res
    integer :: nc
    call run_case( 5, 2, [ ( 0.2_dp, nc = 1, 5 ) ], 0.05_dp, 10, res, nc )
    call verdict( 'constant objective terminates', feasible( res%w ) .and. res%passes <= 10, nf )
  end subroutine

  ! NaN must not propagate into the vector, and must not hang the search.
  subroutine case_nan_objective( nf )
    integer, intent(inout) :: nf
    type(cd_result_ty) :: res
    integer :: nc
    call run_case( 5, 3, [ ( 0.2_dp, nc = 1, 5 ) ], 0.05_dp, 5, res, nc )
    call verdict( 'NaN objective leaves a finite vector', feasible( res%w ), nf )
  end subroutine

  ! Every evaluation reports failure: the search must stop, not spin.
  subroutine case_all_fail( nf )
    integer, intent(inout) :: nf
    type(cd_result_ty) :: res
    integer :: nc
    call run_case( 5, 4, [ ( 0.2_dp, nc = 1, 5 ) ], 0.05_dp, 5, res, nc )
    call verdict( 'all evaluations failing terminates', feasible( res%w ), nf )
  end subroutine

  ! A step larger than the entire simplex: every candidate clamps to a bound.
  subroutine case_huge_delta( nf )
    integer, intent(inout) :: nf
    type(cd_result_ty) :: res
    integer :: nc
    call run_case( 5, 1, [ ( 0.2_dp, nc = 1, 5 ) ], 2.0_dp, 5, res, nc )
    call verdict( 'delta larger than the simplex', feasible( res%w ), nf )
  end subroutine

  ! A step far below the objective's resolution: must still terminate.
  subroutine case_tiny_delta( nf )
    integer, intent(inout) :: nf
    type(cd_result_ty) :: res
    integer :: nc
    call run_case( 5, 1, [ ( 0.2_dp, nc = 1, 5 ) ], 1.0e-9_dp, 3, res, nc )
    call verdict( 'delta far below resolution', feasible( res%w ), nf )
  end subroutine

  subroutine case_one_pass( nf )
    integer, intent(inout) :: nf
    type(cd_result_ty) :: res
    integer :: nc
    call run_case( 5, 1, [ ( 0.2_dp, nc = 1, 5 ) ], 0.05_dp, 1, res, nc )
    call verdict( 'max_pass = 1', feasible( res%w ), nf )
  end subroutine

  ! Input that does not sum to 1 must be normalized, not rejected silently.
  subroutine case_unnormalized_input( nf )
    integer, intent(inout) :: nf
    type(cd_result_ty) :: res
    integer :: nc
    call run_case( 4, 1, [ 3.0_dp, 1.0_dp, 1.0_dp, 1.0_dp ], 0.05_dp, 5, res, nc )
    call verdict( 'unnormalized input is normalized', feasible( res%w ), nf )
  end subroutine

  ! The unconstrained optimum lies far outside the simplex: the search must
  ! stop at the boundary rather than chase it.
  subroutine case_optimum_outside( nf )
    integer, intent(inout) :: nf
    type(cd_result_ty) :: res
    integer :: nc
    call run_case( 5, 5, [ ( 0.2_dp, nc = 1, 5 ) ], 0.05_dp, 10, res, nc )
    call verdict( 'optimum outside the simplex', feasible( res%w ), nf )
  end subroutine

end program test_edge
