#!/usr/bin/env bash
# ADVISORY ONLY -- never a ctest gate, and deliberately not wired into one.
#
# Wall-time smoke check for the hybrid decomposition: hold ranks x threads constant and shift the
# split from all-OpenMP to all-MPI. A correct build gives roughly comparable times; an order of
# magnitude apart means one axis is not parallelising (e.g. every rank integrating every particle).
#
# It is advisory because a timing gate on a shared machine is a flaky gate: this host has 96 cores and
# no exclusivity, so a single number here means little. Read the trend, not the values. Correctness
# across rank counts is gated for real by tests/mpi/, which is the thing that must never regress.
#
# Usage: mpi_scaling_smoke.sh [case-reldir] [total-slots]
set -uo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
CASE="${1:-breakup/khrt-stress}"        # heaviest case in the suite: ~701 children
SLOTS="${2:-8}"
MPIEXEC="${MPIEXEC:-mpiexec}"

command -v "$MPIEXEC" >/dev/null || { echo "no mpiexec on PATH; skipping"; exit 0; }
if ! ldd "$ROOT/bin/IGLOO" 2>/dev/null | grep -q libmpi; then
    echo "bin/IGLOO is not an MPI build (bin/IGLOO is one link target shared by every build"
    echo "tree -- rebuild the MPI tree first). Skipping."
    exit 0
fi

WORK=$(mktemp -d)
trap 'rm -rf "$WORK"' EXIT
cp -rL "$ROOT/tests/$CASE/INPUT" "$WORK/" 2>/dev/null
for f in "$ROOT/tests/$CASE"/*.ini; do cp -L "$f" "$WORK/"; done
[ -d "$ROOT/tests/$CASE/MESH" ] && cp -rL "$ROOT/tests/$CASE/MESH" "$WORK/"
cd "$WORK" || exit 1

printf "%s, %s slots (advisory; wall time on a shared host)\n" "$CASE" "$SLOTS"
printf "%6s %8s %10s\n" "ranks" "threads" "seconds"
n=1
while [ "$n" -le "$SLOTS" ]; do
    t=$(( SLOTS / n ))
    [ "$t" -ge 1 ] || break
    rm -rf OUTPUT; mkdir -p OUTPUT
    #> Nanoseconds as an integer: `bc` is not installed on sprop2, so no floating-point shell math.
    start=$(date +%s%N)
    OMP_NUM_THREADS=$t KMP_STACKSIZE=100M \
        "$MPIEXEC" -n "$n" "$ROOT/bin/IGLOO" >run.txt 2>err.txt
    rc=$?
    end=$(date +%s%N)
    if [ $rc -ne 0 ]; then
        printf "%6s %8s %10s\n" "$n" "$t" "rc=$rc"
    else
        printf "%6s %8s %7d.%01d\n" "$n" "$t" $(( (end-start)/1000000000 )) \
                                            $(( ((end-start)/100000000) % 10 ))
    fi
    n=$(( n * 2 ))
done
