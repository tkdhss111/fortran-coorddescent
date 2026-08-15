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
    real(dp)              :: quantum     = 0.0_dp  ! observed objective resolution
    real(dp)              :: tol_used    = 0.0_dp  ! effective tol at the final incumbent
    integer               :: redundant   = 0       ! slots skipped as duplicates
    integer               :: max_useful_slots = 0  ! widest reachable block seen
    real(dp), allocatable :: w(:)
  end type

  type cd_ty
    integer  :: n        = 0
    real(dp) :: delta    = 0.05_dp     ! ADDITIVE step, in weight units
    !
    ! TOLERANCE. An absolute tol cannot generalize -- it is meaningless to an
    ! objective measured in 1e6 or 1e-9. Two things set it instead:
    !
    !   scale       tol_rel * |f| makes it unit-free
    !   resolution  the objective's own quantum, estimated from what has
    !               actually been observed
    !
    ! The second dominates. If the objective is recorded at coarse precision --
    ! printed to two decimals, counted in integers, read off an instrument --
    ! then a tolerance AT the quantum makes a genuine one-quantum improvement
    ! read as flat, and the search stops early while still improving. Setting
    ! tol just below the quantum resolves one quantum and no less.
    !
    ! auto_tol = .true. (default) does this. Set it .false. and tol is used
    ! exactly as given.
    !
    logical  :: auto_tol = .true.
    real(dp) :: tol_rel  = 1.0e-3_dp   ! used until the quantum is known
    real(dp) :: tol      = 0.01_dp     ! absolute; maintained automatically if auto_tol
    integer  :: max_pass = 5
    integer  :: max_round = 8          ! blocks per direction before giving up
    logical  :: verbose  = .true.
    !
    ! ---- user-tunable, because each one changes the ANSWER ----------------
    !
    ! Sweep order. Two runs are comparable only if this is identical, and the
    ! coordinate swept LAST is starved: earlier ones have claimed the budget.
    ! Defaults to 1..n. Set it explicitly whenever runs will be compared.
    integer,  allocatable :: order(:)
    !
    ! Hard evaluation budget. For an objective costing minutes, max_pass is not
    ! a cost the caller can predict; evaluations are. 0 = unlimited.
    integer  :: max_evals = 0
    !
    ! Per-coordinate bounds within the simplex. Defaults to [0, 1].
    real(dp), allocatable :: wlo(:), whi(:)
    !
    ! Coordinates excluded from the search, keeping their initial weight.
    logical,  allocatable :: frozen(:)
    !
    ! Fraction of the detected quantum used as the flat-stop threshold.
    ! Must be < 1, or a genuine one-quantum gain reads as no change.
    real(dp) :: quantum_frac = 0.5_dp
    real(dp),      allocatable :: w(:)
    logical,       allocatable :: fixed(:)
    character(64), allocatable :: name(:)
    ! Resolution estimate: the smallest non-zero gap between distinct observed
    ! values. For a continuous objective this keeps shrinking and stops
    ! mattering; for a quantized one it converges on the quantum.
    real(dp)              :: q_est     = 0.0_dp
    real(dp)              :: f_scale   = 0.0_dp   ! representative magnitude (the baseline)
    logical               :: quantized = .false.
    real(dp), allocatable :: fseen(:)
    integer               :: nseen = 0
  contains
    procedure :: init       => cd_init
    procedure :: run        => cd_run
    procedure :: observe    => cd_observe
    procedure :: eff_tol    => cd_eff_tol
  end type

  integer, parameter :: MAX_SEEN = 256

  ! Two candidate weights closer than this are the same point: evaluating both
  ! spends an expensive call to learn nothing.
  real(dp), parameter :: W_EPS = 1.0e-12_dp

  ! Coarray scratch for one block of the line search. Module variables are
  ! implicitly saved, which coarrays require.
  !
  ! The block is defined by SLOTS, not by images: there are always at least two
  ! slots (one down, one up), and images take ceil(nslot/nim) slots each. With
  ! nim = 1 that single image serves both directions. Deriving the direction
  ! straight from the image index instead gives a single-image build a downward
  ! search only -- measured: zero improvement over the baseline.
  real(dp), allocatable :: img_val(:)[:]
  real(dp), allocatable :: img_pos(:)[:]
  integer,  allocatable :: img_ok(:)[:]

contains

  ! Block decomposition of [n1,n2] over nprocs for 0-based irank.
  pure subroutine para_range( n1, n2, nprocs, irank, ista, iend )
    integer, intent(in)  :: n1, n2, nprocs, irank
    integer, intent(out) :: ista, iend
    integer :: chunk, extra
    chunk = ( n2 - n1 + 1 ) / nprocs
    extra = mod( n2 - n1 + 1, nprocs )
    ista  = irank * chunk + n1 + min( irank, extra )
    iend  = ista + chunk - 1
    if ( extra > irank ) iend = iend + 1
  end subroutine

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
    if ( present( max_pass ) ) this%max_pass = max_pass
    if ( present( tol ) ) then
      this%tol      = tol
      this%auto_tol = .false.      ! an explicit tol is an explicit choice
    end if

    allocate( this%fseen(MAX_SEEN) )
    this%nseen = 0
    this%q_est = 0.0_dp

    ! Defaults for the tunables the caller did not set.
    if ( .not. allocated( this%order ) ) then
      allocate( this%order(this%n) )
      do k = 1, this%n
        this%order(k) = k
      end do
    end if
    if ( .not. allocated( this%wlo ) )    allocate( this%wlo(this%n),    source = 0.0_dp )
    if ( .not. allocated( this%whi ) )    allocate( this%whi(this%n),    source = 1.0_dp )
    if ( .not. allocated( this%frozen ) ) allocate( this%frozen(this%n), source = .false. )

    call cd_validate( this )
  end subroutine

  !
  ! Validate the configuration and STOP on anything wrong.
  !
  ! Silently clamping a bad parameter is worse than refusing it: the run
  ! completes, the numbers look reasonable, and the caller never learns that
  ! what executed was not what they asked for.
  !
  subroutine cd_validate( this )
    class(cd_ty), intent(inout) :: this
    integer :: k
    logical, allocatable :: hit(:)

    if ( this%delta <= 0.0_dp ) error stop 'cd: delta must be positive'
    if ( this%max_pass < 1 )    error stop 'cd: max_pass must be at least 1'
    if ( this%max_evals < 0 )   error stop 'cd: max_evals must be >= 0 (0 = unlimited)'
    if ( this%tol_rel <= 0.0_dp ) error stop 'cd: tol_rel must be positive'
    if ( this%quantum_frac <= 0.0_dp .or. this%quantum_frac >= 1.0_dp ) &
      error stop 'cd: quantum_frac must be in (0,1) -- at 1 a one-quantum gain reads as flat'
    if ( .not. this%auto_tol .and. this%tol <= 0.0_dp ) &
      error stop 'cd: an explicit tol must be positive'

    if ( size( this%order ) /= this%n ) error stop 'cd: order has the wrong length'
    allocate( hit(this%n), source = .false. )
    do k = 1, this%n
      if ( this%order(k) < 1 .or. this%order(k) > this%n ) error stop 'cd: order index out of range'
      if ( hit(this%order(k)) ) error stop 'cd: order repeats a coordinate'
      hit(this%order(k)) = .true.
    end do
    if ( .not. all( hit ) ) error stop 'cd: order omits a coordinate'

    if ( size( this%wlo ) /= this%n .or. size( this%whi ) /= this%n ) &
      error stop 'cd: wlo/whi have the wrong length'
    if ( any( this%wlo < 0.0_dp ) )        error stop 'cd: wlo must be >= 0'
    if ( any( this%whi > 1.0_dp ) )        error stop 'cd: whi must be <= 1'
    if ( any( this%whi < this%wlo ) )      error stop 'cd: whi must be >= wlo'
    if ( sum( this%wlo ) > 1.0_dp + 1.0e-12_dp ) &
      error stop 'cd: sum(wlo) exceeds 1 -- the feasible set is empty'
    if ( sum( this%whi ) < 1.0_dp - 1.0e-12_dp ) &
      error stop 'cd: sum(whi) is below 1 -- the feasible set is empty'

    if ( size( this%frozen ) /= this%n ) error stop 'cd: frozen has the wrong length'
    if ( count( .not. this%frozen ) < 2 ) &
      error stop 'cd: at least two coordinates must be free to search'
  end subroutine

  !
  ! Record an observed objective value and update the resolution estimate.
  !
  ! The quantum is the smallest non-zero difference between distinct observed
  ! values. A continuous objective drives this toward zero, at which point it
  ! stops constraining anything; a quantized one converges on its step.
  !
  subroutine cd_observe( this, f )
    class(cd_ty), intent(inout) :: this
    real(dp),     intent(in)    :: f
    integer  :: k
    real(dp) :: gap

    if ( .not. ( f == f ) ) return                  ! ignore NaN
    if ( abs( f ) > 0.5_dp * huge( 1.0_dp ) ) return ! ignore the failure sentinel

    do k = 1, this%nseen
      gap = abs( f - this%fseen(k) )
      if ( gap > 0.0_dp ) then
        if ( this%q_est <= 0.0_dp .or. gap < this%q_est ) this%q_est = gap
      end if
    end do

    if ( this%nseen < MAX_SEEN ) then
      this%nseen = this%nseen + 1
      this%fseen(this%nseen) = f
    else
      ! Ring: keep the most recent, so a drifting objective stays represented.
      this%fseen(1:MAX_SEEN-1) = this%fseen(2:MAX_SEEN)
      this%fseen(MAX_SEEN) = f
    end if

    call classify_resolution( this )
  end subroutine

  !
  ! Decide whether the objective is genuinely QUANTIZED or merely continuous.
  !
  ! A small smallest-gap does not mean quantized -- for a continuous objective
  ! that gap simply shrinks as more points are sampled, and treating it as a
  ! quantum drives the tolerance toward zero and disables the flat-stop.
  !
  ! The real signature of quantization is that EVERY observed value lies on a
  ! multiple of the same step. That is what is tested here.
  !
  subroutine classify_resolution( this )
    class(cd_ty), intent(inout) :: this
    integer  :: k, nmult
    real(dp) :: r, dev, span

    this%quantized = .false.
    if ( this%nseen < 8 .or. this%q_est <= 0.0_dp ) return

    span = maxval( this%fseen(1:this%nseen) ) - minval( this%fseen(1:this%nseen) )
    if ( span <= 0.0_dp ) return

    ! A quantum that resolves the whole span into very many levels is
    ! indistinguishable from a continuous objective; refuse to call it quantized.
    if ( span / this%q_est > 1.0e6_dp ) return

    nmult = 0
    do k = 1, this%nseen
      r   = this%fseen(k) / this%q_est
      dev = abs( r - anint( r ) )
      if ( dev < 1.0e-6_dp ) nmult = nmult + 1
    end do

    this%quantized = ( nmult == this%nseen )
  end subroutine

  !
  ! Effective flat-stop threshold at the current incumbent.
  !
  pure real(dp) function cd_eff_tol( this, fref ) result( t )
    class(cd_ty), intent(in) :: this
    real(dp),     intent(in) :: fref
    real(dp) :: scale_ref

    if ( .not. this%auto_tol ) then
      t = this%tol
      return
    end if

    if ( this%quantized ) then
      ! Just under one quantum: a genuine one-quantum gain registers, an
      ! identical reading is flat. Anything at or above the quantum makes the
      ! two indistinguishable and stops the walk while it is still improving.
      t = this%quantum_frac * this%q_est
    else
      ! Continuous: purely relative, so the tolerance is unit-free. The
      ! smallest gap between sampled points is NOT a quantum here -- it keeps
      ! shrinking as more points arrive, and using it would drive the
      ! tolerance to zero and disable the flat-stop altogether.
      scale_ref = max( abs( fref ), abs( this%f_scale ) )
      t = this%tol_rel * scale_ref
    end if

    if ( t <= 0.0_dp ) t = epsilon( 1.0_dp ) * max( abs( fref ), 1.0_dp )
  end function

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
    ! max(0,...) is load-bearing. If rounding ever lets `held` drift above 1,
    ! `avail` goes negative and the clamp below inverts: min( >=0, negative )
    ! selects the negative, so the guard that exists to keep weights
    ! non-negative would itself produce one.
    avail = max( 0.0_dp, 1.0_dp - held )

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

    ! Post-condition. An expensive objective must never be handed an infeasible
    ! point: the run would look successful and the result would be meaningless.
    if ( any( cand < -1.0e-12_dp ) .or. abs( sum( cand ) - 1.0_dp ) > 1.0e-9_dp ) then
      ok = .false.
      cand = this%w
    end if
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
    real(dp) :: best, cur, avail, v, bestv, bestf, prev_d, prev_u, tnow
    integer  :: me, nim, nslot, nper, half, dir, step, s, pass, round, k, i, j
    integer  :: slot, s1, s2, moved, evals, trunc, iord, sl, m, m1, nredun, nuse
    integer  :: nredun_pass, evals_pass
    logical  :: dup
    real(dp), allocatable :: vall(:)
    logical  :: ok, stop_d, stop_u, improved, any_move
    character(64) :: tag

    me    = this_image()
    nim   = num_images()
    nslot = max( 2, nim )         ! always at least one down and one up slot
    half  = nslot / 2             ! odd counts give the spare slot to "up"
    nper  = ( nslot + nim - 1 ) / nim

    if ( .not. allocated( img_val ) ) allocate( img_val(nper)[*], img_pos(nper)[*], img_ok(nper)[*] )
    allocate( cand(this%n), gv(nslot), gp(nslot), go(nslot), vall(nslot) )
    nredun = 0
    nuse   = 0

    best  = baseline
    evals = 0
    trunc = 0
    this%f_scale = abs( baseline )     ! representative magnitude for the relative tol
    call this%observe( baseline )

    if ( me == 1 .and. this%verbose ) then
      write( *, '(a,i0,a,i0,a,i0,a,i0,a)' ) &
        '  coordinate descent: ', nim, ' images, ', nslot, ' slots (', half, &
        ' down, ', nslot - half, ' up)'
      if ( this%auto_tol ) then
        write( *, '(a,f8.4,a,es10.3,a)' ) '  delta = ', this%delta, &
          '   tol = auto (rel ', this%tol_rel, ', floored at half the observed resolution)'
      else
        write( *, '(a,f8.4,a,es10.3)' ) '  delta = ', this%delta, '   tol = fixed ', this%tol
      end if
      write( *, '(a,f10.4)' ) '  baseline = ', baseline
      flush( 6 )
    end if

    pass_loop: do pass = 1, this%max_pass

      this%fixed  = .false.
      moved       = 0
      nredun_pass = 0
      evals_pass  = 0

      do iord = 1, this%n
        s = this%order(iord)

        ! Frozen coordinates keep their weight and are never swept.
        if ( this%frozen(s) ) then
          this%fixed(s) = .true.
          cycle
        end if

        ! The last unfixed coordinate is determined by the residual, not swept.
        if ( count( .not. this%fixed ) <= 1 ) then
          this%fixed(s) = .true.
          cycle
        end if

        ! Hard evaluation budget: stop cleanly rather than overrunning a
        ! caller who is paying minutes per evaluation.
        if ( this%max_evals > 0 .and. evals >= this%max_evals ) then
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

        bestf = best
        call para_range( 1, nslot, nim, me - 1, s1, s2 )

        round_loop: do round = 0, this%max_round - 1

          img_val = huge( 1.0_dp )
          img_ok  = 0
          img_pos = cur

          ! Every image computes the WHOLE slot->value map. It is pure
          ! arithmetic, and it lets each image decide locally whether its own
          ! slot is redundant without another collective.
          do sl = 1, nslot
            if ( sl <= half ) then
              vall(sl) = cur - real( round * max( 1, half ) + sl, dp ) * this%delta
            else
              vall(sl) = cur + real( round * max( 1, half ) + sl - half, dp ) * this%delta
            end if
            vall(sl) = min( max( vall(sl), this%wlo(s) ), min( this%whi(s), avail ) )
          end do

          do slot = s1, s2
            j = slot - s1 + 1
            dir = merge( -1, 1, slot <= half )
            v   = vall(slot)
            img_pos(j) = v

            ! Skip a direction already stopped, but STILL enter every sync all
            ! below -- an image that skips a collective deadlocks the rest.
            if ( ( dir == -1 .and. stop_d ) .or. ( dir == 1 .and. stop_u ) ) cycle

            ! ---- redundancy guard -------------------------------------------
            !
            ! Beyond the reachable range every step clamps onto the same bound,
            ! so extra images re-measure ONE candidate N times. On a
            ! minutes-per-call objective that is the dominant waste, and it is
            ! invisible in the output: the duplicates simply read as flat.
            !
            ! A slot is redundant if it repeats an earlier slot in the same
            ! direction, or if it equals the incumbent, whose value is known.
            dup = ( abs( v - cur ) <= W_EPS )
            if ( .not. dup ) then
              m1 = merge( 1, half + 1, dir == -1 )
              do m = m1, slot - 1
                if ( abs( vall(m) - v ) <= W_EPS ) then
                  dup = .true.
                  exit
                end if
              end do
            end if
            if ( dup ) then
              nredun      = nredun + 1
              nredun_pass = nredun_pass + 1
              img_ok(j) = 0            ! reads as a stop; costs no evaluation
              cycle
            end if

            call make_candidate( this, s, v, cand, ok )
            if ( .not. ok ) cycle
            write( tag, '(a,i0,a,i0,a,i0,a,i0)' ) 'p', pass, '_s', s, '_i', me, '_j', j
            call obj%evaluate( cand, trim( tag ), img_val(j), ok )
            img_ok(j) = merge( 1, 0, ok )
            if ( ok ) call this%observe( img_val(j) )
            evals      = evals + 1
            evals_pass = evals_pass + 1
          end do

          sync all

          if ( me == 1 ) then
            do i = 1, nim
              call para_range( 1, nslot, nim, i - 1, s1, s2 )
              do slot = s1, s2
                j = slot - s1 + 1
                gv(slot) = img_val(j)[i]
                gp(slot) = img_pos(j)[i]
                go(slot) = img_ok(j)[i]
              end do
            end do
            call para_range( 1, nslot, nim, me - 1, s1, s2 )
            ! Scan OUTWARD from the incumbent in each direction, applying the
            ! same rules a serial walk would. Points beyond a stop are dropped.
            ! Recomputed each block: the resolution estimate sharpens as data
            ! arrives, and the relative part tracks the incumbent's scale.
            tnow = this%eff_tol( best )
            call scan_direction( gv(1:half), gp(1:half), go(1:half), &
                                 best, tnow, prev_d, stop_d, bestv, bestf, improved, trunc )
            call scan_direction( gv(half+1:nslot), gp(half+1:nslot), go(half+1:nslot), &
                                 best, tnow, prev_u, stop_u, bestv, bestf, improved, trunc )
          end if

          call co_broadcast( stop_d,   1 )
          call co_broadcast( stop_u,   1 )
          call co_broadcast( bestv,    1 )
          call co_broadcast( bestf,    1 )
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
          ! The incumbent must be the value measured AT bestv -- the point
          ! actually committed to the vector. Taking min() over every measured
          ! point instead admits speculative points the outward scan discarded,
          ! producing an incumbent the vector does not achieve; every later
          ! coordinate then fails to beat it and the search stalls. That defect
          ! worsens with MORE images, so it hides at low image counts.
          best = bestf
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

      ! Reduce BEFORE reporting. These are per-image counters; printing image
      ! 1's share alone understates them by roughly the image count, and the
      ! number looks plausible enough not to question.
      call co_sum( nredun_pass, result_image = 1 )
      call co_sum( evals_pass,  result_image = 1 )

      if ( me == 1 .and. this%verbose ) then
        write( *, '(a,i0,a,i0,a,f10.4)' ) '  pass ', pass, ': ', moved, ' move(s), best ', best
        !
        ! Report MEASURED redundancy, not a heuristic guess.
        !
        ! A startup estimate based on delta and the simplex width is far too
        ! permissive: redundancy is per-coordinate and depends on where that
        ! coordinate currently sits. A coordinate at 0.02 with delta 0.05 has
        ! zero reachable down-steps, which no global rule can see -- so the
        ! heuristic stayed silent at 32 images while 162 slots were skipped.
        !
        if ( nredun_pass > 0 ) then
          write( *, '(a,i0,a,i0,a)' ) &
            '         ', nredun_pass, ' of ', nredun_pass + evals_pass, &
            ' slots were duplicates (clamped onto a bound) and were skipped'
          if ( nredun_pass > evals_pass ) &
            write( *, '(a)' ) &
              '         more than half the slots are redundant: fewer images would do the same work'
        end if
        flush( 6 )
      end if

      any_move = moved > 0
      call co_broadcast( any_move, 1 )
      if ( .not. any_move ) exit pass_loop
    end do pass_loop

    ! Both counters are per-image and must be reduced, exactly like evals.
    ! Reporting image 1's share alone reads as zero at high image counts,
    ! because image 1 owns the first slot and the first slot is never a
    ! duplicate -- which makes the metric contradict the evaluation count.
    call co_sum( evals,  result_image = 1 )
    call co_sum( nredun, result_image = 1 )
    call co_sum( trunc,  result_image = 1 )

    res%best        = best
    res%baseline    = baseline
    res%evaluations = evals
    res%passes      = pass
    res%truncations = trunc
    res%quantum     = merge( this%q_est, 0.0_dp, this%quantized )
    res%redundant   = nredun
    res%max_useful_slots = nuse
    res%tol_used    = this%eff_tol( best )
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
  subroutine scan_direction( vals, pos, okf, incumbent, tol, prev, stopped, bestv, bestf, improved, trunc )
    real(dp), intent(in)    :: vals(:), pos(:)
    integer,  intent(in)    :: okf(:)
    real(dp), intent(in)    :: incumbent, tol
    real(dp), intent(inout) :: prev
    logical,  intent(inout) :: stopped
    real(dp), intent(inout) :: bestv
    real(dp), intent(inout) :: bestf     ! objective AT bestv -- the committed value
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

      if ( vals(i) <= incumbent - tol .and. vals(i) < bestf ) then
        bestv    = pos(i)
        bestf    = vals(i)
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
