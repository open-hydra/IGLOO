#!/bin/bash
#SBATCH --job-name=igloo-final-vie-plait-dense
#SBATCH --partition=free
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=80
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=02:00:00
#
# The scaling curve of vie-plait-dense (compute-bound, gas order 2) on one exclusive node, IGLOO_REPS repetitions
# (default 5): the 1x1 anchor, every split of 4, 8, 16 and 32 slots, the 80-slot splits, the A/A twin.
# Submit from the harness root:
#   sbatch --output=$IGLOO_PERF_ROOT/logs/%x-%j.out --export=ALL,IGLOO_PERF_ROOT=$IGLOO_PERF_ROOT \
#     scripts/igloo-final-vie-plait-dense.slurm
set -u
W=${IGLOO_PERF_ROOT:?}/final
source "$IGLOO_PERF_ROOT/scripts/_igloo_common.sh"
igloo_setup || exit 1
bash "$ROOT/scripts/make-workloads.sh" vie-plait-dense || { echo "IGLOO-FINAL STOPPED: the workload did not build"; exit 1; }
REPS=${IGLOO_REPS:-5}

igloo_matrix vie-plait-dense "$REPS" check

echo "IGLOO-FINAL DONE"
