#!/bin/bash
#SBATCH --job-name=igloo-final
#SBATCH --partition=free
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=80
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=16:00:00
#
# The scaling matrix on one exclusive node: per workload the 1x1 anchor, every split of 4, 8, 16 and 32 slots and
# the 80-slot splits, IGLOO_REPS repetitions (default 5); the A/A twin; khrt-stress with OUTPUT/ on tmpfs and with
# the scatter cloud off. Every point is referred to the 1x1 anchor of its workload in this job (analyze.py).
# Submit from the harness root, after the smoke:
#   sbatch --output=$IGLOO_PERF_ROOT/logs/%x-%j.out --export=ALL,IGLOO_PERF_ROOT=$IGLOO_PERF_ROOT scripts/igloo-final.slurm
set -u
W=${IGLOO_PERF_ROOT:?}/final
source "$IGLOO_PERF_ROOT/scripts/_igloo_common.sh"
igloo_setup || exit 1
bash "$ROOT/scripts/make-workloads.sh" || { echo "IGLOO-FINAL STOPPED: the workloads did not build"; exit 1; }
REPS=${IGLOO_REPS:-5}

for wl in khrt-stress vie-plait-dense box-4fam; do
  chk=check; [ "$wl" = box-4fam ] && chk=finite
  echo "=== $wl: the 1x1 anchor ==="
  igloo_run "$wl-1x1" "$wl" 1 1 1 "$(igloo_spread 1)" "$REPS" "$chk"

  echo "=== $wl: every split of 4, 8, 16 and 32 slots, the slots spread over the sockets ==="
  for N in 4 8 16 32; do
    for R in 1 2 4 8 16 32; do
      [ "$R" -le "$N" ] || continue
      Tn=$(( N / R ))
      igloo_run "$wl-${R}x${Tn}" "$wl" "$R" "$Tn" "$R" "$(igloo_spread "$N")" "$REPS"
    done
  done

  echo "=== $wl: the 80-slot splits, Intel MPI's omp domains ==="
  for RT in 1x80 4x20 8x10 16x5 20x4 40x2 80x1; do
    # 80 ranks of the box: 80 x 2.09 GiB at the finalize peak, beside the node's memory
    [ "$wl" = box-4fam ] && [ "$RT" = 80x1 ] && continue
    R=${RT%x*}; Tn=${RT#*x}
    igloo_run "$wl-$RT" "$wl" "$R" "$Tn" "$R" NONE "$REPS"
  done
done

echo "=== the A/A twin: a second snapshot of the same binary ==="
for wl in khrt-stress vie-plait-dense; do
  IGLOO_AA=1 igloo_run "$wl-1x1-aa" "$wl" 1 1 1 "$(igloo_spread 1)" "$REPS"
done
for wl in khrt-stress vie-plait-dense box-4fam; do
  IGLOO_AA=1 igloo_run "$wl-4x20-aa" "$wl" 4 20 4 NONE "$REPS"
done

echo "=== khrt-stress: OUTPUT/ on tmpfs ==="
IGLOO_OUTPUT=shm igloo_run khrt-stress-1x80-shm khrt-stress  1 80  1 NONE "$REPS"
IGLOO_OUTPUT=shm igloo_run khrt-stress-4x20-shm khrt-stress  4 20  4 NONE "$REPS"
IGLOO_OUTPUT=shm igloo_run khrt-stress-80x1-shm khrt-stress 80  1 80 NONE "$REPS"

echo "=== khrt-stress: the scatter cloud off, source and Eulerian output (the accumulators' atomics) ==="
IGLOO_INI_SET="out-scatter=off out-file=ALL" \
  igloo_run khrt-stress-1x1-noscat  khrt-stress 1  1 1 "$(igloo_spread 1)"  "$REPS"
IGLOO_INI_SET="out-scatter=off out-file=ALL" \
  igloo_run khrt-stress-1x32-noscat khrt-stress 1 32 1 "$(igloo_spread 32)" "$REPS"

echo "IGLOO-FINAL DONE"
