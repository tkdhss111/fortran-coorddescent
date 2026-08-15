# fortran-coorddescent

Simplex-constrained coordinate descent with a **coarray-parallel line search** —
for expensive black-box objectives (minutes per evaluation) subject to
`sum(w) = 1, w >= 0`.

## What is parallelized, and why it is exact

The line search along a single coordinate is where a serial descent spends
everything, so that is what gets distributed. With N images the lower half
search **downward** and the upper half **upward**, each image taking a different
number of steps:

```
N = 4    image 1 -> cur - 1*delta      image 3 -> cur + 1*delta
         image 2 -> cur - 2*delta      image 4 -> cur + 2*delta

N = 8    images 1-4 -> cur - {1,2,3,4}*delta
         images 5-8 -> cur + {1,2,3,4}*delta
```

The result is **identical to a serial walk**, not an approximation. The stop
rules depend on the *ordering* of the points, not on the order in which they
were measured, so after `sync all` image 1 scans outward from the incumbent and
applies exactly the serial rules, discarding anything past a stop. Speculating
costs evaluations; it never costs correctness.

If a direction needs more steps than N/2, another block runs. Images whose
direction has already stopped still enter every `sync all` — skipping a
collective on one image deadlocks the rest.

An odd image count gives the spare image to the upward direction.

## Stop rules

Both are required; either ends that direction.

```
flat        |f(k) - f(k-1)| <= tol
worsening   f(k) > incumbent  .and.  f(k) >= f(k-1)
```

The second is not optional. Without it a walk that has passed the optimum keeps
stepping away from it, and on a unimodal coordinate every further point is
guaranteed worse. Measured waste before it was added: **40–44 % of all
evaluations**.

## Truncations are a result, not bookkeeping

`cd_result_ty%truncations` counts walks ended by the *flat* rule while still
improving — capped by the tolerance rather than by an optimum. A later pass
would continue them.

Report this number. Reading a tolerance-capped coordinate as "converged" is how
a search silently understates a weight, and the reader cannot tell the two apart
from the final vector alone.

## Withdraw-and-redistribute

Coordinate `s` takes value `v`; the mass held by **fixed** coordinates is
withdrawn, and what remains after `s` takes its share is redistributed among the
still-unfixed coordinates in proportion to their current values. The constraint
therefore holds *exactly* at every candidate, rather than being restored by a
projection afterwards.

Where the remaining pool is all zero the residual is shared equally, so the map
stays well defined at the corners of the simplex — which is exactly where a
sparse optimum lives.

**Known limitation.** By the time the sweep reaches the last few coordinates,
earlier ones have claimed most of the budget, so `R = 1 - sum(fixed)` is small
and the late coordinates get almost no search. Multiple passes mitigate it: each
pass unfixes everything, letting starved coordinates compete from the new
vector. Sweep order is therefore not neutral to the outcome — fix it explicitly
if you intend to compare two runs.

## API

```fortran
use coord_descent_mo, only: cd_ty, cd_result_ty, objective_ty

type(cd_ty)        :: cd
type(cd_result_ty) :: res

call cd%init ( w0, names, delta = 0.05_dp, tol = 0.01_dp, max_pass = 10 )
call cd%run  ( obj, baseline, res )

! res%w  res%best  res%evaluations  res%passes  res%truncations
```

The objective is supplied by the caller as a type extending `objective_ty`:

```fortran
subroutine evaluate ( this, w, tag, val, ok )
```

It is a type rather than a procedure pointer because an expensive objective
usually carries state — paths, a run counter, a cache — and every image calls it
concurrently. `tag` is unique per candidate, so an implementation that shells out
to an external process can give each evaluation its own working directory.

`ok = .false.` marks a failed evaluation; the search treats it as a stop rather
than as a good point. A crashed run must never look like an improvement.

## Choosing delta

`delta` is **additive**, in weight units, deliberately. A multiplicative grid
makes the step proportional to the current value, so a coordinate that starts
small can never take a large step — from a uniform start every coordinate reads
"flat" and nothing moves, which looks exactly like convergence.

Pick it against the resolution of your objective. Near an optimum the response is
second order:

```
f(w) ~ f* + 0.5 * k * (w - w*)^2      so a step h moves f by 0.5*k*h^2
```

so a step must satisfy `0.5*k*h^2 > tol` to be visible at all. If your objective
is recorded at coarse precision, check that before choosing anything — a
tolerance larger than the response makes every optimizer report "no move"
regardless of the weights.

## Building

```sh
cmake -S . -B build -G Ninja -DCMAKE_Fortran_COMPILER=ifx \
      -DCOARRAY_MODE=shared -DCOARRAY_IMAGES=4
cmake --build build

source /opt/intel/oneapi/setvars.sh      # required to RUN, not only to build
./build/test_coord_descent
```

`-coarray=shared` launches images through MPI Hydra, so oneAPI must be sourced
before running. Profile with `-coarray=single`: a shared build is already running
N images, so profiling it measures the wrong thing.

## Validation

`test/test_coord_descent.f90` minimizes a quadratic with a **known** optimum on
the simplex, deliberately sparse so the corner behaviour is exercised:

```
f(w) = sum_k a_k * (w_k - t_k)^2      minimized at w = t
```

Result at 8 coordinates, 4 images, delta 0.05:

```
converged at pass 9, 0 moves
f                 0.456750 -> 0.003252
max |w - t|       0.0341        (< one delta)
sum(w)            1.000000      exactly
sparse targets    recovered to exactly 0.0000
evaluations       326
tolerance-capped  25
```

Give it enough passes. At `max_pass = 6` the same run stops mid-descent — still
making 3 moves per pass — and strands a coordinate at 0.1 whose true value is 0.
The final vector looks plausible either way; only the move count at the last
pass distinguishes "converged" from "ran out of passes".

## Licence

MIT.
