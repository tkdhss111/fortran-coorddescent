#!/usr/bin/env bash
#
# Run the correctness test across image counts.
#
# A coarray program that passes at one image count proves nothing about the
# others. The original version of this module passed at 2, 3 and 4 images and
# was BROKEN at 1, 6 and 8 -- a single-image build searched downward only, and
# at higher counts an over-optimistic incumbent stalled the search. Neither was
# visible from the 4-image run.
#
# Run this, not a single configuration.
#
set -uo pipefail
BIN=${1:-./build/test_coord_descent}
COUNTS=${2:-"1 2 3 4 5 6 8 12 16"}
fail=0

printf "  images  rc   final_f      max|w-t|   evals   converged\n"
for n in $COUNTS; do
  out=$(FOR_COARRAY_NUM_IMAGES=$n timeout 900 "$BIN" 2>&1)
  rc=$?
  final=$(echo "$out" | grep "final f"     | awk '{print $3}')
  err=$(echo   "$out" | grep "max |w - t|" | awk '{print $NF}')
  ev=$(echo    "$out" | grep "evaluations" | awk '{print $2}')
  conv=$(echo  "$out" | grep -oE "pass [0-9]+: 0 move" | head -1)
  printf "  %-7s %-4s %-12s %-10s %-7s %s\n" \
    "$n" "$rc" "${final:-NONE}" "${err:-NONE}" "${ev:-NONE}" "${conv:-NOT CONVERGED}"
  [ "$rc" -ne 0 ] && fail=$((fail+1))
  [ -z "$conv" ] && { echo "    ^ did not converge: raise max_pass"; fail=$((fail+1)); }
done

echo
if [ "$fail" -eq 0 ]; then
  echo "  ALL IMAGE COUNTS PASSED"
else
  echo "  $fail CONFIGURATION(S) FAILED"
  exit 1
fi
