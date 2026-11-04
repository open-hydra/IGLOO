#!/bin/bash
#SBATCH --job-name=igloo-numa
#SBATCH --partition=free
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=80
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=02:00:00
#
# Page placement per workload: one rank of 80 and of 40 threads over the sockets, pages where the first touch left
# them (plain) and interleaved (inter), beside the four-rank layout; IGLOO_REPS repetitions (default 2). No anchor:
# the question is the ratio inside this job. Submit from the harness root:
#   sbatch --output=$IGLOO_PERF_ROOT/logs/%x-%j.out --export=ALL,IGLOO_PERF_ROOT=$IGLOO_PERF_ROOT scripts/igloo-numa.slurm
set -u
W=${IGLOO_PERF_ROOT:?}/numa
source "$IGLOO_PERF_ROOT/scripts/_igloo_common.sh"
igloo_setup || exit 1
bash "$ROOT/scripts/make-workloads.sh" || { echo "IGLOO-NUMA STOPPED: the workloads did not build"; exit 1; }
REPS=${IGLOO_REPS:-2}

for wl in khrt-stress vie-plait-dense box-4fam; do
  for Tn in 80 40; do
    echo "=== $wl 1 x $Tn: plain, then interleaved ==="
    igloo_run "$wl-1x$Tn-plain" "$wl" 1 "$Tn" 1 "$(igloo_spread "$Tn")" "$REPS"
    IGLOO_WRAP="numactl --interleave=all" \
      igloo_run "$wl-1x$Tn-inter" "$wl" 1 "$Tn" 1 "$(igloo_spread "$Tn")" "$REPS"
  done
  echo "=== $wl: the four-rank layout on the same node ==="
  igloo_run "$wl-4x20" "$wl" 4 20 4 NONE "$REPS"
done

echo "IGLOO-NUMA DONE"
