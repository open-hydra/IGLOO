#!/bin/bash
#SBATCH --job-name=igloo-final-khrt-stress
#SBATCH --partition=free
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=80
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=02:00:00
#
# The scaling curve of khrt-stress (output-bound) on one exclusive node, IGLOO_REPS repetitions (default 5): the
# 1x1 anchor, every split of 4, 8, 16 and 32 slots, the 80-slot splits, the A/A twin; then OUTPUT/ on tmpfs, and
# the scatter cloud off with source and Eulerian output (the accumulators' atomics). Submit from the harness root:
#   sbatch --output=$IGLOO_PERF_ROOT/logs/%x-%j.out --export=ALL,IGLOO_PERF_ROOT=$IGLOO_PERF_ROOT \
#     scripts/igloo-final-khrt-stress.slurm
set -u
W=${IGLOO_PERF_ROOT:?}/final
source "$IGLOO_PERF_ROOT/scripts/_igloo_common.sh"
igloo_setup || exit 1
bash "$ROOT/scripts/make-workloads.sh" khrt-stress || { echo "IGLOO-FINAL STOPPED: the workload did not build"; exit 1; }
REPS=${IGLOO_REPS:-5}

igloo_matrix khrt-stress "$REPS" check

echo "=== khrt-stress: OUTPUT/ on tmpfs ==="
IGLOO_OUTPUT=shm igloo_run khrt-stress-1x80-shm khrt-stress  1 80  1 NONE "$REPS"
IGLOO_OUTPUT=shm igloo_run khrt-stress-4x20-shm khrt-stress  4 20  4 NONE "$REPS"
IGLOO_OUTPUT=shm igloo_run khrt-stress-80x1-shm khrt-stress 80  1 80 NONE "$REPS"

echo "=== khrt-stress: the scatter cloud off, source and Eulerian output ==="
IGLOO_INI_SET="out-scatter=off out-file=ALL" \
  igloo_run khrt-stress-1x1-noscat  khrt-stress 1  1 1 "$(igloo_spread 1)"  "$REPS"
IGLOO_INI_SET="out-scatter=off out-file=ALL" \
  igloo_run khrt-stress-1x32-noscat khrt-stress 1 32 1 "$(igloo_spread 32)" "$REPS"

echo "IGLOO-FINAL DONE"
