#!/usr/bin/env bash
# The time wrapper of igloo-io.slurm: runs its arguments and prints one ITIMES line on stderr with their user and
# system CPU time, tagged with the MPI rank (GNU time when installed, else the shell's `times` of its children).
#   IGLOO_WRAP="bash io-times.sh [env NAME=value ...]"     (igloo_run puts the binary after the wrapper)
r=${PMI_RANK:-${PMIX_RANK:-${OMPI_COMM_WORLD_RANK:-${SLURM_PROCID:-x}}}}
if [ -x /usr/bin/time ]; then
  exec /usr/bin/time -f "ITIMES rank=$r user_s=%U sys_s=%S wall_s=%e maxrss_kb=%M fs_in=%I fs_out=%O vcsw=%w ivcsw=%c rc=%x" "$@"
fi
"$@"
rc=$?
f=$(mktemp "${TMPDIR:-/tmp}/io-times.XXXXXX") || exit $rc
times > "$f"
awk -v r="$r" -v rc="$rc" '
  function sec(t,  m) { m = index(t, "m"); return substr(t, 1, m - 1) * 60 + substr(t, m + 1, length(t) - m - 1) }
  NR == 2 { printf "ITIMES rank=%s user_s=%.2f sys_s=%.2f rc=%s source=times\n", r, sec($1), sec($2), rc }' "$f" >&2
rm -f "$f"
exit $rc
