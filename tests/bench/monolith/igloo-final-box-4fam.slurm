#!/bin/bash
#SBATCH --job-name=igloo-final-box-4fam
#SBATCH --partition=free
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=80
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=20:00:00
#
# The scaling curve of box-4fam (the reduction volume and the memory per rank) on one exclusive node, IGLOO_REPS
# repetitions (default 3: every run pays rank 0's serial read of the ASCII gas file and its write of the grid
# output): the 1x1 anchor, every split of 4, 8, 16 and 32 slots, the 80-slot splits but 80 x 1, the A/A twin at
# 4x20. Submit from the harness root:
#   sbatch --output=$IGLOO_PERF_ROOT/logs/%x-%j.out --export=ALL,IGLOO_PERF_ROOT=$IGLOO_PERF_ROOT \
#     scripts/igloo-final-box-4fam.slurm
set -u
W=${IGLOO_PERF_ROOT:?}/final
source "$IGLOO_PERF_ROOT/scripts/_igloo_common.sh"
igloo_setup || exit 1
bash "$ROOT/scripts/make-workloads.sh" box-4fam || { echo "IGLOO-FINAL STOPPED: the workload did not build"; exit 1; }
REPS=${IGLOO_REPS:-3}

igloo_matrix box-4fam "$REPS" finite

echo "IGLOO-FINAL DONE"
