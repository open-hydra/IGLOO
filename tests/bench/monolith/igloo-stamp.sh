#!/usr/bin/env bash
# One rank of a timed launch: a stamp when the rank starts, the command, a stamp when it ends, one file per rank in
# $IGLOO_STAMPS. igloo_run launches `mpirun ... bash igloo-stamp.sh <command ...>` and reads the stamps: the time
# from its own clock before mpirun to the first start (launch_s), from the first start to the last end (life_s),
# and from the last end to mpirun's return (teardown_s).
r=${PMI_RANK:-${PMIX_RANK:-${OMPI_COMM_WORLD_RANK:-${SLURM_PROCID:-x}}}}
f=${IGLOO_STAMPS:?IGLOO_STAMPS: the stamp directory}/rank.$r
echo "start $(date +%s.%N)" >> "$f"
"$@"
rc=$?
echo "end $(date +%s.%N) $rc" >> "$f"
exit $rc
