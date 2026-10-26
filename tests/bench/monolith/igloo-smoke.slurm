#!/bin/bash
#SBATCH --job-name=igloo-smoke
#SBATCH --partition=free
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=80
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=01:00:00
#
# Pre-flight of the IGLOO campaign: the binary runs on a compute node, the pinning takes, the workloads build, and
# each workload runs once at 4 ranks x 20 threads (one rank per socket) with its checks. Submit from the harness root:
#   sbatch --output=$IGLOO_PERF_ROOT/logs/%x-%j.out --export=ALL,IGLOO_PERF_ROOT=$IGLOO_PERF_ROOT scripts/igloo-smoke.slurm
set -u
W=${IGLOO_PERF_ROOT:?}/smoke
source "$IGLOO_PERF_ROOT/scripts/_igloo_common.sh"
igloo_setup || exit 1
echo "perf_event_paranoid $(cat /proc/sys/kernel/perf_event_paranoid 2>/dev/null || echo unknown)"
echo "PERF $(command -v perf || echo none) NUMACTL $(command -v numactl || echo none)"

t0=$(date +%s)
bash "$ROOT/scripts/make-workloads.sh"; rc=$?
echo "IMAKE rc=$rc s=$(( $(date +%s) - t0 ))"
[ "$rc" -eq 0 ] || { echo "IGLOO-SMOKE STOPPED: the workloads did not build"; exit 1; }

igloo_run khrt-stress-4x20     khrt-stress     4 20 4 NONE 1 check
igloo_run vie-plait-dense-4x20 vie-plait-dense 4 20 4 NONE 1 check
igloo_run box-4fam-4x20        box-4fam        4 20 4 NONE 1 finite

echo "IGLOO-SMOKE DONE"
