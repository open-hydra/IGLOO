#!/bin/bash
#SBATCH --job-name=igloo-omp
#SBATCH --partition=free
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=80
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=05:00:00
#
# Thread placement per workload: one rank inside one socket (NUMA-free), the same thread counts spread over the
# sockets (the cost of remote pages), one rank of 40 and 80 threads plain and with the pages interleaved, and the
# hybrid controls 4x20 and 16x5; IGLOO_REPS repetitions (default 2). Submit from the harness root:
#   sbatch --output=$IGLOO_PERF_ROOT/logs/%x-%j.out --export=ALL,IGLOO_PERF_ROOT=$IGLOO_PERF_ROOT scripts/igloo-omp.slurm
set -u
W=${IGLOO_PERF_ROOT:?}/omp
source "$IGLOO_PERF_ROOT/scripts/_igloo_common.sh"
igloo_setup || exit 1
bash "$ROOT/scripts/make-workloads.sh" || { echo "IGLOO-OMP STOPPED: the workloads did not build"; exit 1; }
REPS=${IGLOO_REPS:-2}

for wl in khrt-stress vie-plait-dense box-4fam; do
  echo "=== $wl A: one rank inside socket 0 (the 1x1 point is the anchor) ==="
  for Tn in 1 2 4 5 10 20; do
    tag="$wl-1x$Tn-s0"; [ "$Tn" -eq 1 ] && tag="$wl-1x1"
    igloo_run "$tag" "$wl" 1 "$Tn" 1 "$(igloo_sock0 "$Tn")" "$REPS"
  done
  echo "=== $wl B: the same thread counts spread over the sockets ==="
  for Tn in 4 8 20 40; do
    igloo_run "$wl-1x$Tn-sp" "$wl" 1 "$Tn" 1 "$(igloo_spread "$Tn")" "$REPS"
  done
  echo "=== $wl C: one rank of 40 and 80 threads, pages where the first touch left them, then interleaved ==="
  for Tn in 40 80; do
    igloo_run "$wl-1x$Tn-plain" "$wl" 1 "$Tn" 1 "$(igloo_spread "$Tn")" "$REPS"
    IGLOO_WRAP="numactl --interleave=all" \
      igloo_run "$wl-1x$Tn-inter" "$wl" 1 "$Tn" 1 "$(igloo_spread "$Tn")" "$REPS"
  done
  echo "=== $wl D: the hybrid controls ==="
  igloo_run "$wl-4x20" "$wl"  4 20  4 NONE "$REPS"
  igloo_run "$wl-16x5" "$wl" 16  5 16 NONE "$REPS"
done

echo "IGLOO-OMP DONE"
