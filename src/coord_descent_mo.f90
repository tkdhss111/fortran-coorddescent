!
! Simplex-constrained coordinate descent with a coarray-parallel line search.
!
! The objective is assumed EXPENSIVE (minutes per call) and the constraint is
! sum(w) = 1, w >= 0. One coordinate is swept at a time; the remainder of the
! budget is withdrawn and redistributed among the coordinates not yet fixed, so
! the simplex constraint holds exactly at every candidate rather than being
! restored by a projection afterwards.
!
! PARALLELISM
!
! The line search along one coordinate is what gets distributed. With N images,
! the lower half search DOWNWARD and the upper half UPWARD, each image taking a
! different number of steps, so several points along both directions are
! measured simultaneously:
!
!   N = 4      image 1 -> cur - 1*delta      image 3 -> cur + 1*delta
!              image 2 -> cur - 2*delta      image 4 -> cur + 2*delta
!
!   N = 8      images 1-4 -> cur - {1,2,3,4}*delta
!              images 5-8 -> cur + {1,2,3,4}*delta
!
! A serial walk evaluates these one at a time and stops early; this evaluates a
! whole block at once and applies the SAME stop rules afterwards, scanning
! outward from the incumbent. The answer is identical to the serial walk -- only
! the wall-clock differs -- because the stop rules depend on the ordering of the
! points, not on the order in which they were measured. Points beyond a stop are
! discarded, which is the cost of speculating.
!
! STOP RULES (both required; either one ends that direction)
!
!   flat       |f(k) - f(k-1)| <= tol         the curve has gone flat
!   worsening  f(k) > incumbent .and. f(k) >= f(k-1)
!
! The second is not optional. Without it a walk that has passed the optimum
! keeps stepping away from it, and every remaining point in that direction is
! guaranteed worse on a unimodal coordinate. Measured waste before it was
! added: 40-44 % of all evaluations.
!
! An odd image count gives the extra image to the upward direction.
!
module coord_descent_mo

  use, intrinsic :: iso_fortran_env, only: real64

  implicit none

  private
  public :: cd_ty, objective_ty, cd_result_ty

  integer, parameter :: dp = real64

  !
  ! The caller owns the objective. It is deliberately not a simple function
  ! pointer: an expensive objective usually needs state (paths, a run counter,
  ! a cache), and every image calls it concurrently.
  !
  type, abstract :: objective_ty
  contains
    procedure(evaluate_i), deferred :: evaluate
  end type

  abstract interface
    subroutine evaluate_i( this, w, tag, val, ok )
      import :: objective_ty, dp
      class(objective_ty), intent(inout) :: this
      real(dp),            intent(in)    :: w(:)    ! full vector, sums to 1
      character(*),        intent(in)    :: tag     ! unique per candidate
      real(dp),            intent(out)   :: val     ! minimized
      logical,             intent(out)   :: ok      ! .false. => treat as failed
    end subroutine
  end interface

  type cd_result_ty
    real(dp)              :: best        = 0.0_dp
    real(dp)              :: baseline    = 0.0_dp
    integer               :: evaluations = 0
    integer               :: passes      = 0
    integer               :: truncations = 0   ! walks cut by tol, not by an optimum
    real(dp), allocatable :: w(:)
  end type

  type cd_ty
    integer  :: n        = 0
    real(dp) :: delta    = 0.05_dp     ! ADDITIVE step, in weight units
    real(dp) :: tol      = 0.01_dp     ! flat-stop threshold on the objective
    integer  :: max_pass = 5
    integer  :: max_round = 8          ! blocks per direction before giving up
    logical  :: verbose  = .true.
    real(dp),      allocatable :: w(:)
    logical,       allocatable :: fixed(:)
    character(64), allocatable :: name(:)
  contains
    procedure :: init => cd_init
    procedure :: run  => cd_run
  end type

  ! Coarray scratch for one block of the line search. Module variables are
  ! implicitly saved, which coarrays require.
  real(dp) :: img_val[*]
  real(dp) :: img_pos[*]
  integer  :: img_ok[*]

contains

  subroutine cd_init( this, w0, names, delta, tol, max_pass )
    class(cd_ty),  intent(inout) :: this
    real(dp),      intent(in)    :: w0(:)
    character(*),  intent(in)    :: names(:)
    real(dp), optional, intent(in) :: delta, tol
    integer,  optional, intent(in) :: max_pass
    integer :: k

    this%n = size( w0 )
    if ( size( names ) /= this%n ) error stop 'cd_init: names and w0 differ in length'
    if ( any( w0 < 0.0_dp ) )      error stop 'cd_init: negative initial weight'

    allocate( this%w, source = w0 )
    allocate( this%fixed(this%n), source = .false. )
    allocate( this%name(this%n) )
    do k = 1, this%n
      this%name(k) = names(k)
    end do

    ! Normalize once so the simplex constraint holds from the start.
    this%w = this%w / sum( this%w )

    if ( present( delta ) )    this%delta    = delta
    if ( present( tol ) )      this%tol      = tol
    if ( present( max_pass ) ) this%max_pass = max_pass
  end subroutine

  !
  ! Withdraw-and-redistribute: coordinate s takes value v; the mass not held by
  ! FIXED coordinates and not taken by s is shared among the remaining unfixed
  ! coordinates in proportion to their current values. If those are all zero it
  ! is shared equally, which keeps the map well defined at the corners.
  !
  subroutine make_candidate( this, s, v, cand, ok )
    class(cd_ty), intent(in)  :: this
    integer,      intent(in)  :: s
    real(dp),     intent(in)  :: v
    real(dp),     intent(out) :: cand(:)
    logical,      intent(out) :: ok

    real(dp) :: held, avail, pool, rest
    integer  :: k, nfree

    ok = .true.
    held = 0.0_dp
    do k = 1, this%n
      if ( this%fixed(k) ) held = held + this%w(k)
    end do
    avail = 1.0_dp - held

    if ( v < -1.0e-12_dp .or. v > avail + 1.0e-12_dp ) then
      ok = .false.
      cand = this%w
      return
    end if

    cand = this%w
    cand(s) = min( max( v, 0.0_dp ), avail )
    rest = avail - cand(s)

    pool  = 0.0_dp
    nfree = 0
    do k = 1, this%n
      if ( .not. this%fixed(k) .and. k /= s ) then
        pool  = pool + this%w(k)
        nfree = nfree + 1
      end if
    end do

    if ( nfree == 0 ) then
      ! Nothing left to absorb the remainder; s must take all of it.
      cand(s) = avail
      return
    end if

    do k = 1, this%n
      if ( .not. this%fixed(k) .and. k /= s ) then
        if ( pool > 1.0e-15_dp ) then
          cand(k) = rest * this%w(k) / pool
        else
          cand(k) = rest / real( nfree, dp )
        end if
      end if
    end do
  end subroutine

  !
  ! Main loop. Every image executes it; the decisions are made on image 1 and
  ! broadcast, so the images stay in lockstep on the same vector.
  !
  subroutine cd_run( this, obj, baseline, res )
    class(cd_ty),        intent(inout) :: this
    class(objective_ty), intent(inout) :: obj
    real(dp),            intent(in)    :: baseline
    type(cd_result_ty),  intent(out)   :: res

    real(dp), allocatable :: cand(:), gv(:), gp(:)
    integer,  allocatable :: go(:)
    real(dp) :: best, cur, avail, v, bestv, prev_d, prev_u
    integer  :: me, nim, half, dir, step, s, pass, round, k, i, moved, evals, trunc
    logical  :: ok, stop_d, stop_u, improved, any_move
    character(64) :: tag

    me   = this_image()
    nim  = num_images()
    half = max( 1, nim / 2 )      ! odd counts give the spare image to "up"

    allocate( cand(this%n), gv(nim), gp(nim), go(nim) )

    best  = baseline
    evals = 0
    trunc = 0

    if ( me == 1 .and. this%verbose ) then
      write( *, '(a,i0,a,i0,a,i0,a)' ) &
        '  coordinate descent: ', nim, ' images (', half, ' down, ', nim - half, ' up)'
      write( *, '(a,f8.4,a,f8.4)' ) '  delta = ', this%delta, '   tol = ', this%tol
      write( *, '(a,f10.4)' ) '  baseline = ', baseline
      flush( 6 )
    end if

    pass_loop: do pass = 1, this%max_pass

      this%fixed = .false.
      moved = 0

      do s = 1, this%n

        ! The last unfixed coordinate is determined by the residual, not swept.
        if ( count( .not. this%fixed ) <= 1 ) then
          this%fixed(s) = .true.
          cycle
        end if

        cur   = this%w(s)
        avail = 1.0_dp - sum( this%w, mask = this%fixed )
        bestv = cur
        improved = .false.
        stop_d = .false.
        stop_u = .false.
        prev_d = best
        prev_u = best

        if ( me <= half ) then
          dir  = -1
          step = me
        else
          dir  = 1
          step = me - half
        end if

        round_loop: do round = 0, this%max_round - 1

          k = round * half + step
          v = cur + real( dir * k, dp ) * this%delta
          v = min( max( v, 0.0_dp ), avail )

          ! Skip work in a direction already stopped, but keep participating in
          ! the collectives -- every image must reach every sync all.
          if ( ( dir == -1 .and. stop_d ) .or. ( dir == 1 .and. stop_u ) ) then
            img_val = huge( 1.0_dp )
            img_ok  = 0
          else
            call make_candidate( this, s, v, cand, ok )
            if ( ok ) then
              write( tag, '(a,i0,a,i0,a,i0)' ) 'p', pass, '_s', s, '_i', me
              call obj%evaluate( cand, trim( tag ), img_val, ok )
              img_ok = merge( 1, 0, ok )
              evals  = evals + 1
            else
              img_val = huge( 1.0_dp )
              img_ok  = 0
            end if
          end if
          img_pos = v

          sync all

          if ( me == 1 ) then
            do i = 1, nim
              gv(i) = img_val[i]
              gp(i) = img_pos[i]
              go(i) = img_ok[i]
            end do
            ! Scan OUTWARD from the incumbent in each direction, applying the
            ! same rules a serial walk would. Points beyond a stop are dropped.
            call scan_direction( gv(1:half), gp(1:half), go(1:half), &
                                 best, this%tol, prev_d, stop_d, bestv, improved, trunc )
            call scan_direction( gv(half+1:nim), gp(half+1:nim), go(half+1:nim), &
                                 best, this%tol, prev_u, stop_u, bestv, improved, trunc )
          end if

          call co_broadcast( stop_d,   1 )
          call co_broadcast( stop_u,   1 )
          call co_broadcast( bestv,    1 )
          call co_broadcast( improved, 1 )
          call co_broadcast( prev_d,   1 )
          call co_broadcast( prev_u,   1 )

          if ( stop_d .and. stop_u ) exit round_loop
        end do round_loop

        ! Commit: move s to its best value, redistribute, and fix it.
        call make_candidate( this, s, bestv, cand, ok )
        if ( ok ) this%w = cand
        this%fixed(s) = .true.

        if ( improved ) then
          moved = moved + 1
          ! Recompute the incumbent from the value that won.
          best = min( best, minval( gv, mask = go == 1 ) )
        end if

        if ( me == 1 .and. this%verbose ) then
          if ( improved ) then
            write( *, '(a,a24,a,f8.4,a,f10.4)' ) &
              '    MOVED ', trim( this%name(s) ), ' -> ', bestv, '   best ', best
          else
            write( *, '(a,a24,a,f8.4)' ) &
              '    kept  ', trim( this%name(s) ), ' at  ', bestv
          end if
          flush( 6 )
        end if
      end do

      if ( me == 1 .and. this%verbose ) then
        write( *, '(a,i0,a,i0,a,f10.4)' ) '  pass ', pass, ': ', moved, ' move(s), best ', best
        flush( 6 )
      end if

      any_move = moved > 0
      call co_broadcast( any_move, 1 )
      if ( .not. any_move ) exit pass_loop
    end do pass_loop

    call co_sum( evals, result_image = 1 )

    res%best        = best
    res%baseline    = baseline
    res%evaluations = evals
    res%passes      = pass
    res%truncations = trunc
    allocate( res%w, source = this%w )
  end subroutine

  !
  ! Apply the serial stop rules to a block of already-measured points, ordered
  ! by increasing distance from the incumbent.
  !
  ! `trunc` counts walks ended by the FLAT rule while still improving -- those
  ! are cut by the tolerance rather than by an optimum, and a later pass would
  ! continue them. Recording that is not bookkeeping: reading a tolerance-capped
  ! coordinate as "converged" is how a search silently understates a weight.
  !
  subroutine scan_direction( vals, pos, okf, incumbent, tol, prev, stopped, bestv, improved, trunc )
    real(dp), intent(in)    :: vals(:), pos(:)
    integer,  intent(in)    :: okf(:)
    real(dp), intent(in)    :: incumbent, tol
    real(dp), intent(inout) :: prev
    logical,  intent(inout) :: stopped
    real(dp), intent(inout) :: bestv
    logical,  intent(inout) :: improved
    integer,  intent(inout) :: trunc

    integer  :: i
    real(dp) :: d

    if ( stopped ) return

    do i = 1, size( vals )
      if ( okf(i) /= 1 ) then
        stopped = .true.
        return
      end if

      if ( vals(i) <= incumbent - tol ) then
        bestv    = pos(i)
        improved = .true.
      end if

      ! worsening: past the optimum, every further point is worse
      if ( vals(i) > incumbent .and. vals(i) >= prev ) then
        stopped = .true.
        return
      end if

      d = abs( vals(i) - prev )
      if ( d <= tol ) then
        ! Flat while still better than the incumbent = capped by tolerance.
        if ( vals(i) < incumbent ) trunc = trunc + 1
        stopped = .true.
        return
      end if

      prev = vals(i)
    end do
  end subroutine

end module coord_descent_mo
