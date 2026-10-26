#!/usr/bin/env bash
# Bring a job back from the cluster: its .out file and its per-run logs (the job's WORKDIR/logs/<jobid>/), into
# <dest>/<jobid>/ (dest default $IGLOO_PERF_LOCAL). A build job has no WORKDIR: only its .out comes back.
#   IGLOO_PERF_REMOTE=<ssh host> IGLOO_PERF_ROOT=<harness root there> fetch-results.sh <jobid> [<dest>]
set -eu
JOB=${1:?job id}
REMOTE=${IGLOO_PERF_REMOTE:?IGLOO_PERF_REMOTE: the ssh host of the cluster}
ROOT=${IGLOO_PERF_ROOT:?IGLOO_PERF_ROOT: the harness root on the cluster}
DEST=${2:-${IGLOO_PERF_LOCAL:?IGLOO_PERF_LOCAL: the local directory of the fetched jobs}}/$JOB
mkdir -p "$DEST/logs"

rsync -a "$REMOTE:$ROOT/logs/*-$JOB.out" "$DEST/"
out=$(ls "$DEST"/*-"$JOB".out | head -1)
W=$(sed -n 's/^WORKDIR //p' "$out" | head -1)
if [ -n "$W" ]; then
  rsync -a "$REMOTE:$W/logs/$JOB/" "$DEST/logs/"
fi
echo "fetched $(ls "$DEST/logs" | wc -l) run logs and $(basename "$out") into $DEST"
