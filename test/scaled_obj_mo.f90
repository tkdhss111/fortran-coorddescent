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
