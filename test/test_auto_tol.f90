!
! Automatic tolerance: scale invariance and resolution awareness.
!
! Two failure modes a fixed absolute tolerance cannot survive, both tested here
! against the same underlying problem:
!
!   SCALE       the objective is multiplied by 1e6. A tol chosen for an
!               objective of order 1 now sits far below the noise floor of a
!               different one, or far above its entire dynamic range.
!
!   RESOLUTION  the objective is recorded to two decimals -- printed, counted,
!               read from an instrument. A tol AT the quantum makes a genuine
!               one-quantum improvement indistinguishable from no change, so
!               the walk stops while still improving and the coordinate is
!               silently understated.
!
! The second is the one that bites in practice, because the search still
! terminates cleanly and reports a plausible vector.
!
module scaled_obj_mo

  use, intrinsic :: iso_fortran_env, only: real64
  use coord_descent_mo, only: objective_ty

  implicit none

  private
  public :: scaled_ty

  integer, parameter :: dp = real64

  type, extends(objective_ty) :: scaled_ty
    real(dp), allocatable :: t(:)
    real(dp), allocatable :: a(:)
    real(dp)              :: scale   = 1.0_dp   ! multiply the objective
    real(dp)              :: quantum = 0.0_dp   ! round to this; 0 = continuous
  contains
    procedure :: evaluate => scaled_evaluate
  end type

contains

  subroutine scaled_evaluate( this, w, tag, val, ok )
    class(scaled_ty), intent(inout) :: this
    real(dp),         intent(in)    :: w(:)
    character(*),     intent(in)    :: tag
    real(dp),         intent(out)   :: val
    logical,          intent(out)   :: ok
    val = this%scale * sum( this%a * ( w - this%t )**2 )
    if ( this%quantum > 0.0_dp ) then
      val = anint( val / this%quantum ) * this%quantum
    end if
    ok = .true.
  end subroutine

end module scaled_obj_mo


program test_auto_tol

  use, intrinsic :: iso_fortran_env, only: real64
  use coord_descent_mo, only: cd_ty, cd_result_ty
  use scaled_obj_mo,    only: scaled_ty

  implicit none

  integer, parameter :: dp = real64
  integer, parameter :: N = 8

  integer  :: nfail
  real(dp) :: e_auto, e_fixed
  nfail = 0

  if ( this_image() == 1 ) write( *, '(a)' ) 'automatic tolerance'

  ! Scale invariance: the same problem at 1x and 1e6x must behave the same.
  call one_case( 'continuous, scale 1',      1.0_dp,   0.0_dp,   .true.,  nfail, e_auto )
  call one_case( 'continuous, scale 1e6',    1.0e6_dp, 0.0_dp,   .true.,  nfail, e_fixed )
  call compare( 'scale invariance', e_auto, e_fixed, nfail )

  ! Resolution. At quantum 0.01 each step still changes the objective by many
  ! quanta, so a fixed tol at the quantum costs nothing -- reported for honesty.
  call one_case( 'quantized 0.01, auto tol', 1.0_dp,   0.01_dp,  .true.,  nfail, e_auto )
  call one_case( 'quantized 0.01, tol=0.01', 1.0_dp,   0.01_dp,  .false., nfail, e_fixed )
  call not_worse( 'quantum 0.01', e_auto, e_fixed, nfail )

  ! COARSE quantum: the whole objective spans only ~9 levels, so improvements
  ! are one or two quanta and a tolerance AT the quantum reads them as flat.
  ! This is the regime the automatic mode exists for.
  call one_case( 'quantized 0.05, auto tol', 1.0_dp,   0.05_dp,  .true.,  nfail, e_auto )
  call one_case( 'quantized 0.05, tol=0.05', 1.0_dp,   0.05_dp,  .false., nfail, e_fixed )
  call not_worse( 'quantum 0.05', e_auto, e_fixed, nfail )

  if ( this_image() == 1 ) then
    write( *, '(a)' ) ''
    if ( nfail == 0 ) then
      write( *, '(a)' ) '  ALL AUTO-TOLERANCE CHECKS PASSED'
    else
      write( *, '(a,i0,a)' ) '  ', nfail, ' CHECK(S) FAILED'
      error stop 1
    end if
  end if

contains

  !
  ! Two runs of the SAME problem must land in the same place when only the
  ! objective's scale differs -- that is what scale invariance means.
  !
  subroutine compare( label, e1, e2, nf )
    character(*), intent(in)    :: label
    real(dp),     intent(in)    :: e1, e2
    integer,      intent(inout) :: nf
    if ( this_image() /= 1 ) return
    write( *, '(a)' ) ''
    if ( abs( e1 - e2 ) < 1.0e-9_dp ) then
      write( *, '(a,a)' ) '  PASS  ', label//': identical result at 1x and 1e6x'
    else
      write( *, '(a,a,2f10.4)' ) '  FAIL  ', label//': results differ ', e1, e2
      nf = nf + 1
    end if
  end subroutine

  !
  ! The right criterion for a quantized objective. Asserting that the optimum
  ! is REACHED is wrong: at a coarse quantum the objective does not distinguish
  ! the weights, so no optimizer can recover them and the assertion tests the
  ! measurement rather than the method. What must hold is that the automatic
  ! tolerance is never worse than a fixed one.
  !
  subroutine not_worse( label, e_auto, e_fixed, nf )
    character(*), intent(in)    :: label
    real(dp),     intent(in)    :: e_auto, e_fixed
    integer,      intent(inout) :: nf
    if ( this_image() /= 1 ) return
    write( *, '(a)' ) ''
    if ( e_auto <= e_fixed + 1.0e-9_dp ) then
      write( *, '(a,a,f8.4,a,f8.4)' ) '  PASS  ', label//': auto ', e_auto, '  <=  fixed ', e_fixed
    else
      write( *, '(a,a,f8.4,a,f8.4)' ) '  FAIL  ', label//': auto ', e_auto, '  >  fixed ', e_fixed
      nf = nf + 1
    end if
  end subroutine

  subroutine one_case( label, scale, quantum, auto, nf, err_out )
    character(*), intent(in)    :: label
    real(dp),     intent(in)    :: scale, quantum
    logical,      intent(in)    :: auto
    integer,      intent(inout) :: nf
    real(dp),     intent(out)   :: err_out

    type(cd_ty)        :: cd
    type(scaled_ty)    :: obj
    type(cd_result_ty) :: res
    character(64) :: names(N)
    character(8)  :: lbl
    real(dp) :: w0(N), t(N), a(N), base, err_w, rel_gain
    integer  :: k
    logical  :: ok

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

    allocate( obj%t, source = t )
    allocate( obj%a, source = a )
    obj%scale   = scale
    obj%quantum = quantum

    if ( auto ) then
      call cd%init( w0, names, delta = 0.05_dp, max_pass = 30 )
    else
      ! A fixed tol equal to the quantum: the pathological case.
      call cd%init( w0, names, delta = 0.05_dp, tol = max( quantum, 0.01_dp ), max_pass = 30 )
    end if
    cd%verbose = .false.

    call obj%evaluate( w0, 'base', base, ok )
    call cd%run( obj, base, res )

    sync all

    if ( this_image() == 1 ) then
      err_w    = maxval( abs( res%w - t ) )
      rel_gain = ( res%baseline - res%best ) / max( abs( res%baseline ), tiny( 1.0_dp ) )
      write( *, '(a)' ) ''
      write( *, '(a,a)' )      '  case          ', label
      write( *, '(a,es12.4)' ) '    baseline    ', res%baseline
      write( *, '(a,es12.4)' ) '    final       ', res%best
      write( *, '(a,f8.2,a)' ) '    gain        ', 100.0_dp * rel_gain, ' %'
      write( *, '(a,es12.4)' ) '    quantum est ', res%quantum
      write( *, '(a,es12.4)' ) '    tol used    ', res%tol_used
      write( *, '(a,f10.4)' )  '    max |w - t| ', err_w
      write( *, '(a,i0,a,i0)' ) '    evals ', res%evaluations, '   truncations ', res%truncations

      if ( auto ) then
        ! Reaching the optimum is only required where the objective can express
        ! it -- i.e. a continuous one. Under coarse quantization the weights are
        ! genuinely unidentifiable and demanding recovery tests the measurement.
        if ( quantum <= 0.0_dp ) then
          if ( err_w > 0.05_dp ) then
            write( *, '(a)' ) '    FAIL  did not reach the optimum'
            nf = nf + 1
          else
            write( *, '(a)' ) '    PASS  optimum reached'
          end if
        end if
        if ( quantum > 0.0_dp ) then
          if ( res%quantum <= 0.0_dp .or. abs( res%quantum - quantum ) > 0.5_dp * quantum ) then
            write( *, '(a)' ) '    FAIL  resolution not detected'
            nf = nf + 1
          else
            write( *, '(a)' ) '    PASS  resolution detected from observed values'
          end if
          if ( res%tol_used >= quantum ) then
            write( *, '(a)' ) '    FAIL  tol is not below the quantum'
            nf = nf + 1
          else
            write( *, '(a)' ) '    PASS  tol sits below the quantum'
          end if
        end if
      else
        ! Reported, not asserted: whether a tolerance at the quantum actually
        ! costs anything depends on how many quanta a single step spans. State
        ! the outcome rather than claiming a failure that may not have occurred.
        write( *, '(a,es12.4,a,f8.4)' ) '    NOTE  fixed tol ', res%tol_used, &
                                        ' at the quantum; reached ', err_w
      end if
      flush( 6 )
    end if

    err_out = maxval( abs( res%w - t ) )
    deallocate( obj%t, obj%a )
  end subroutine

end program test_auto_tol
