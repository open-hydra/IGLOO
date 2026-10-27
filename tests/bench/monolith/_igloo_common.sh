# Shared functions of the IGLOO timing jobs on the cluster. Sourced, not executed (README.md beside this file).
#
# The protocol: one binary per job, snapshotted into the job's directory; every point of a curve on the node set
# of one job; efficiency against the 1x1 anchor of the same job; every repetition scored, the median and the
# minimum kept with every field of the fastest repetition; resident memory sampled per point; the placement read
# back from the MPI and OpenMP runtimes at launch. The instrument is IGLOO's own [IGLOO-timing] report.
#
# The caller sets W (the job's work directory) before igloo_setup; IGLOO_PERF_ROOT is required. No `set -e` here:
# a failed point must not end the job.

ROOT=${IGLOO_PERF_ROOT:?IGLOO_PERF_ROOT: the harness root on the cluster}
SRC=${IGLOO_PERF_SRC:-$ROOT/src}
CASE=${IGLOO_PERF_CASE:-$ROOT/case}
BIN0=${IGLOO_PERF_BIN:-$ROOT/bin/IGLOO-hyb}
BIN=$BIN0
BIN_AA=$BIN0
STAMP=$ROOT/scripts/igloo-stamp.sh

# Job set-up: limits, pinning defaults, the binary snapshot and its twin, the provenance lines.
igloo_setup () {
  ulimit -s unlimited
  export KMP_STACKSIZE=100M OMP_STACKSIZE=100M
  export I_MPI_PIN_RESPECT_CPUSET=1 I_MPI_PIN_DOMAIN=omp
  export OMP_PROC_BIND=close OMP_PLACES=cores
  [ -x "$BIN0" ] || { echo "IFAIL setup reason=no-binary ($BIN0)"; return 1; }
  [ -f "$STAMP" ] || { echo "IFAIL setup reason=no-stamp-wrapper ($STAMP)"; return 1; }
  LOGS=$W/logs/$SLURM_JOB_ID
  mkdir -p "$LOGS" "$W/work/$SLURM_JOB_ID"
  cd "$W" || return 1
  # The job runs its own copy: a rebuild of $BIN0 during the job cannot splice two binaries into one curve.
  # The copy keeps the basename IGLOO, which the memory sampler and the occupancy check match (ps -C IGLOO).
  mkdir -p "$W/bin.$SLURM_JOB_ID" "$W/bin.$SLURM_JOB_ID-aa"
  cp -f "$BIN0" "$W/bin.$SLURM_JOB_ID/IGLOO" && BIN=$W/bin.$SLURM_JOB_ID/IGLOO
  cp -f "$BIN0" "$W/bin.$SLURM_JOB_ID-aa/IGLOO" && BIN_AA=$W/bin.$SLURM_JOB_ID-aa/IGLOO
  echo "JOB=${SLURM_JOB_NAME:-?} id=$SLURM_JOB_ID nodes=${SLURM_NODELIST:-?} partition=${SLURM_JOB_PARTITION:-?}"
  echo "START $(date -u +%Y-%m-%dT%H:%M:%SZ)"
  echo "WORKDIR $W"
  echo "BINARY md5=$(md5sum < "$BIN" | cut -d' ' -f1) $(stat -c %s "$BIN") bytes (snapshot: $BIN; twin: $BIN_AA)"
  if [ -f "$BIN0.commit" ]; then echo "COMMIT $(tr '\n' ' ' < "$BIN0.commit")"; else echo "COMMIT unknown (no $BIN0.commit)"; fi
  echo "CPU $(lscpu | sed -n 's/^Model name: *//p') | sockets $(lscpu | sed -n 's/^Socket(s): *//p') |" \
       "cores/socket $(lscpu | sed -n 's/^Core(s) per socket: *//p') | threads/core $(lscpu | sed -n 's/^Thread(s) per core: *//p') |" \
       "NUMA $(lscpu | sed -n 's/^NUMA node(s): *//p') | flags-md5 $(lscpu | grep '^Flags' | md5sum | cut -c1-12)"
  free -g | head -2
  echo "LAUNCH RES_OPTIONS=${IGLOO_LAUNCH_RES_OPTIONS:-unset} (IGLOO_LAUNCH_RES_OPTIONS; mpirun only)"
  # The OpenMP wait policy decides how long idle threads spin after a region, and the spin lands in the cycles of
  # the phase that follows: every job states it (the harness sets none of the three).
  echo "WAITPOLICY KMP_BLOCKTIME=${KMP_BLOCKTIME:-unset} OMP_WAIT_POLICY=${OMP_WAIT_POLICY:-unset}" \
       "GOMP_SPINCOUNT=${GOMP_SPINCOUNT:-unset}"
  igloo_workloads
}

# One WORKLOAD line per built workload: its name and the digest of its md5.txt (a new digest = a new baseline).
igloo_workloads () {
  local d
  for d in "$CASE"/*/; do
    [ -f "$d/md5.txt" ] && echo "WORKLOAD $(basename "$d") md5=$(md5sum < "$d/md5.txt" | cut -c1-12)"
  done
  return 0
}

# igloo_domains <core-list> <ranks> <threads>: one hex mask per rank for I_MPI_PIN_DOMAIN. A numactl around
# mpirun does not pin under Slurm (the ranks are spawned through Slurm's bootstrap); Intel MPI pins from these.
igloo_domains () {
  python3 -c '
import sys
cores = [int(x) for x in sys.argv[1].split(",")]
R, T = int(sys.argv[2]), int(sys.argv[3])
if len(cores) != R * T:
    sys.exit("igloo_domains: %d cores for %d ranks x %d threads" % (len(cores), R, T))
out = []
for r in range(R):
    m = 0
    for c in cores[r * T:(r + 1) * T]:
        m |= 1 << c
    out.append("0x%x" % m)
print("[" + ",".join(out) + "]")' "$1" "$2" "$3"
}

# CPU lists from lscpu, one CPU per physical core: igloo_spread N (N cores evenly over the sockets), igloo_sock0 N
# (the first N cores of the first socket), igloo_packed N (the first N cores in order).
igloo_cores () {
  python3 - "$1" "$2" <<'PY'
import collections, subprocess, sys
mode, N = sys.argv[1], int(sys.argv[2])
out = subprocess.run(["lscpu", "-p=CPU,CORE,SOCKET"], stdout=subprocess.PIPE).stdout.decode()
first = {}
for line in out.splitlines():
    if line.startswith("#") or not line.strip():
        continue
    cpu, core, s = (int(x) for x in line.split(",")[:3])
    if (s, core) not in first or cpu < first[(s, core)]:
        first[(s, core)] = cpu
sock = collections.OrderedDict()
for (s, core) in sorted(first):
    sock.setdefault(s, []).append(first[(s, core)])
for s in sock:
    sock[s].sort()
S = len(sock)
if mode == "spread":
    cpus = sorted(sock[list(sock)[i]][0] for i in range(N)) if N < S else \
        sorted(c for s in sock for c in sock[s][:N // S])
elif mode == "sock0":
    cpus = sock[list(sock)[0]][:N]
else:
    cpus = sorted(c for s in sock for c in sock[s])[:N]
if len(cpus) != N:
    sys.exit("igloo_cores: %s %d gives %d cores" % (mode, N, len(cpus)))
print(",".join(str(c) for c in cpus))
PY
}
igloo_spread () { igloo_cores spread "$1"; }
igloo_sock0 () { igloo_cores sock0 "$1"; }
igloo_packed () { igloo_cores packed "$1"; }

# No stale ranks from the previous point: a run that starts beside them measures contention. The process checks
# here and below count this user's IGLOO processes only.
igloo_wait_clear () {
  local i n=0
  for i in $(seq 1 60); do
    n=$(pgrep -u "$(id -u)" -x IGLOO 2>/dev/null | wc -l)
    [ "$n" -eq 0 ] && return 0
    [ "$i" -eq 1 ] && echo "IWARN stale IGLOO ranks=$n -- waiting"
    sleep 5
  done
  echo "IWARN stale IGLOO ranks=$n -- PROCEEDING ANYWAY"
}

# igloo_rmshm <dir>: removes a tmpfs output directory and its job directory once empty (no-op on "").
igloo_rmshm () {
  [ -n "$1" ] || return 0
  rm -rf "$1"; rmdir "${1%/*}" 2>/dev/null
  return 0
}

# Where the point runs, sampled once while it runs: threads named IGLOO and the distinct cores they sit on.
igloo_occupancy () {
  ps -eLo psr,uid,comm --no-headers 2>/dev/null | awk -v tag="$1" -v uid="$(id -u)" '
    $3=="IGLOO" && $2==uid { sol[$1]++; ns++ }
    END { printf "IOCC tag=%s solver_threads=%d solver_cores=%d\n", tag, ns, length(sol) }'
}

# igloo_ini_set <input.ini> key=value ...: each key set in [IGLOO-General] (its old line replaced, else inserted
# after the section header).
igloo_ini_set () {
  python3 - "$@" <<'PY'
import re, sys
path, pairs = sys.argv[1], [a.split("=", 1) for a in sys.argv[2:]]
if any(len(p) != 2 or not p[0] for p in pairs):
    sys.exit("igloo_ini_set: arguments are key=value")
keys = [k.strip() for k, _ in pairs]
lines = open(path).read().split("\n")
out, sec, found = [], None, False
for line in lines:
    s = line.strip()
    if s.startswith("["):
        sec = s
        out.append(line)
        if s == "[IGLOO-General]":
            found = True
            out.extend("%s = %s" % (k.strip(), v.strip()) for k, v in pairs)
        continue
    if sec == "[IGLOO-General]" and any(re.match(re.escape(k) + r"\s*=", s) for k in keys):
        continue
    out.append(line)
if not found:
    sys.exit("igloo_ini_set: no [IGLOO-General] section in %s" % path)
open(path, "w").write("\n".join(out))
PY
}

# igloo_casedir <dir> <workload>: a run directory the workload's own checks can run in. <dir> stands for the
# tests/ directory: every entry of $SRC/tests but the workload's category is a symlink, the category is a real
# directory holding symlinks to its other cases, and the case itself holds INPUT (a symlink to the workload's),
# input.ini and check.py (copies) and an empty OUTPUT/. Sets CASEDIR.
igloo_casedir () {
  local D=$1 wl=$2 cat e
  [ -f "$CASE/$wl/md5.txt" ] || { echo "IFAIL workload=$wl reason=workload-not-built ($CASE/$wl)"; return 1; }
  cat=$(cat "$CASE/$wl/category.txt")
  rm -rf "$D"
  mkdir -p "$D/$cat/$wl/OUTPUT"
  for e in "$SRC"/tests/*; do
    [ "${e##*/}" = "$cat" ] || ln -s "$e" "$D/${e##*/}"
  done
  if [ -d "$SRC/tests/$cat" ]; then
    for e in "$SRC/tests/$cat"/*; do
      [ "${e##*/}" = "$wl" ] || ln -s "$e" "$D/$cat/${e##*/}"
    done
  fi
  ln -s "$CASE/$wl/INPUT" "$D/$cat/$wl/INPUT"
  cp "$CASE/$wl/input.ini" "$D/$cat/$wl/input.ini"
  [ -f "$CASE/$wl/check.py" ] && cp "$CASE/$wl/check.py" "$D/$cat/$wl/check.py"
  CASEDIR=$D/$cat/$wl
  return 0
}

# igloo_parse <tag> <rep> <stdout> <stderr>: line 1 the IAFF line (placement at launch: I_MPI_DEBUG's rank table,
# OMP_DISPLAY_AFFINITY's thread lines), line 2 the report's fields as key=value, with the cycles and the wait
# policy of the launch.
igloo_parse () {
  python3 - "$@" <<'PY'
import os, re, sys
tag, rep, fout, ferr = sys.argv[1:5]
out = open(fout, errors="replace").read()
err = open(ferr, errors="replace").read()
both = out + "\n" + err


def cpus(s):
    got = set()
    for part in s.replace(" ", "").split(","):
        if not part:
            continue
        if "-" in part:
            a, b = part.split("-", 1)
            got.update(range(int(a), int(b) + 1))
        else:
            got.add(int(part))
    return got


def num(s):
    try:
        return float(s)
    except ValueError:
        return float("nan")


# placement: Intel MPI's rank table (I_MPI_DEBUG >= 4) and the OpenMP runtime's per-thread binding
pin = {}
for m in re.finditer(r"MPI startup\(\):\s+(\d+)\s+(\d+)\s+(\S+)\s+\{([^}]*)\}", both):
    pin[int(m.group(1))] = (m.group(3), frozenset(cpus(m.group(4))))
thr = {}
for m in re.finditer(r"pid\s+(\d+)\s+tid\s+(\d+)\s+thread\s+(\d+)\s+bound to OS proc set\s+\{([^}]*)\}", both):
    thr.setdefault(m.group(1), {})[int(m.group(3))] = frozenset(cpus(m.group(4)))
sets = [s for t in thr.values() for s in t.values()]
pin_cpus = set().union(*[c for _, c in pin.values()])
nthr = sorted(len(t) for t in thr.values())
if not nthr:
    threads = "0"
elif nthr[0] == nthr[-1]:
    threads = "%d" % nthr[0]
else:
    threads = "%d-%d" % (nthr[0], nthr[-1])
print("IAFF tag=%s rep=%s ranks=%d threads=%s cores=%d cpus=%d pin_ranks=%d pin_cpus=%d pin_nodes=%s" % (
    tag, rep, len(thr), threads, len(set(sets)), len(set().union(*sets)), len(pin), len(pin_cpus),
    ",".join(sorted({n for n, _ in pin.values()})) or "-"))

f = {}
blocks = out.split("[IGLOO-timing] sweep ")
f["reports"] = len(blocks) - 1
last = blocks[-1] if len(blocks) > 1 else ""
h = re.match(r"(\d+): (\d+) solve\(s\), (\d+) rank\(s\) x (\d+) thread\(s\)", last)
if h:
    f["sweep"], f["solves"], f["hdr_ranks"], f["hdr_threads"] = h.groups()
SETUP = ["setup:input", "setup:gas-read", "setup:blocks", "setup:gas-import", "setup:pinning"]
SWEEP = ["accumulators", "ghost-fill", "solve-other", "parcel-loop", "census", "reduce", "finalize",
         "smoothing", "writeout+merge"]
ph = {}
for m in re.finditer(r"\[IGLOO-timing\] phase (\S+)\s+(\S+) s\s+ranks min\s+(\S+) max\s+(\S+)", last):
    ph[m.group(1)] = tuple(num(m.group(k)) for k in (2, 3, 4))
if ph:
    f["setup_s"] = "%.6f" % sum(ph.get(p, (0.0,))[0] for p in SETUP)
    f["sweep_s"] = "%.6f" % sum(ph.get(p, (0.0,))[0] for p in SWEEP)
    for key, p in (("loop_s", "parcel-loop"), ("reduce_s", "reduce"), ("finalize_s", "finalize"),
                   ("smoothing_s", "smoothing"), ("writeout_s", "writeout+merge")):
        if p in ph:
            f[key] = "%.6f" % ph[p][0]
    for p in SETUP + SWEEP:
        if p in ph:
            f["ph_" + re.sub(r"[^A-Za-z0-9]", "_", p)] = "%.6f/%.6f/%.6f" % ph[p]
rk = re.findall(r"\[IGLOO-timing\] rank\s+(\d+)\s+parcel-loop\s+(\S+) s\s+parcels\s+(\d+)\s+threads\s+(\d+)"
                r"\s+busy\s+(\S+) s\s+thread busy max/mean\s+(\S+)", last)
if rk:
    f["rank_lines"] = len(rk)
    f["parcels"] = sum(int(r[2]) for r in rk)
    f["parcel_loop_max_s"] = "%.6f" % max(num(r[1]) for r in rk)
    f["busy_maxmean_max"] = "%.3f" % max(num(r[5]) for r in rk)
m = re.search(r"\[IGLOO-timing\] imbalance\s+rank parcel-loop max/mean\s+(\S+)", last)
if m:
    f["imbalance"] = m.group(1)
m = re.search(r"\[IGLOO-timing\] sizes per rank \(values\)\s+gas (\d+)\s+geometry (\d+)\s+accumulators (\d+)"
              r"\s+total\s+(\S+) MiB", last)
if m:
    f["sz_gas"], f["sz_geo"], f["sz_acc"], f["sz_mib"] = m.groups()
m = re.search(r"MPI ranks = (\d+)", out)
f["mpi_ranks"] = m.group(1) if m else "1"
n = re.findall(r"number of particles =\s+(\d+)", out)
if n:
    f["particles"] = sum(int(x) for x in n)
n = re.findall(r"number of children =\s+(\d+)", out)
if n:
    f["children"] = sum(int(x) for x in n)
# core cycles (a build with the counter): the state, the totals, every phase; the wait policy they depend on
m = re.search(r"\[IGLOO-timing\] cycles  (counters|unavailable|incomplete)", last)
f["cyc_state"] = m.group(1) if m else "absent"
m = re.search(r"\[IGLOO-timing\] cycles  counters (\d+) of (\d+) thread", last)
if m:
    f["cyc_counters"], f["cyc_threads"] = m.groups()
for m in re.finditer(r"\[IGLOO-timing\] cycles phase (\S+)\s+(\d+)\s+ranks min\s+(\d+) max\s+(\d+)", last):
    f["cy_" + re.sub(r"[^A-Za-z0-9]", "_", m.group(1))] = "%s/%s/%s" % m.group(2, 3, 4)
    if m.group(1) == "parcel-loop":
        f["cyc_loop"] = m.group(2)
m = re.search(r"\[IGLOO-timing\] cycles  sweep (\d+)\s+setup (\d+)\s+aggregate\s+(\S+) GHz", last)
if m:
    f["cyc_sweep"], f["cyc_setup"], f["cyc_ghz"] = m.groups()
for k in ("KMP_BLOCKTIME", "OMP_WAIT_POLICY", "GOMP_SPINCOUNT"):
    f[k.lower()] = os.environ.get(k, "unset").replace(" ", "_") or "unset"
print(" ".join("%s=%s" % kv for kv in f.items()))
PY
}

# igloo_stamps <dir> <t0 ns> <t1 ns>: from igloo-stamp.sh's files, the launch (t0 to the first rank start), the
# ranks' lifetime (first start to last end), the teardown (last end to t1) and the number of ranks stamped.
igloo_stamps () {
  cat "$1"/rank.* 2>/dev/null | awk -v t0="$2" -v t1="$3" '
    $1=="start" { n++; if (s == "" || $2 < s) s = $2 }
    $1=="end"   { if (e == "" || $2 > e) e = $2 }
    END {
      if (n == 0 || e == "") { printf "launch_s= life_s= teardown_s= stamped=%d", n + 0; exit }
      printf "launch_s=%.3f life_s=%.3f teardown_s=%.3f stamped=%d", s - t0 / 1e9, e - s, t1 / 1e9 - e, n }'
}

# igloo_summary <reps file> <reps asked>: the IOK line of a point from its IREP lines -- the median, minimum and
# standard deviation of wall_s, life_s, launch_s, sweep_s, loop_s, reduce_s, cyc_sweep and cyc_loop, then every
# field of the fastest repetition (its rep= field names it).
igloo_summary () {
  python3 - "$1" "$2" <<'PY'
import re, statistics, sys
reps = []
for line in open(sys.argv[1]):
    if line.startswith("IREP "):
        reps.append(dict(re.findall(r"(\w+)=(\S*)", line)))
if not reps:
    sys.exit(0)


def fx(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return float("nan")


best = min(reps, key=lambda d: fx(d.get("wall_s")))
stat = []
for k in ("wall_s", "life_s", "launch_s", "sweep_s", "loop_s", "reduce_s", "cyc_sweep", "cyc_loop"):
    v = [fx(d[k]) for d in reps if d.get(k)]
    if not v:
        continue
    sd = "%.6f" % statistics.stdev(v) if len(v) > 1 else "na"
    stat.append("%s_med=%.6f %s_min=%.6f %s_sd=%s" % (k, statistics.median(v), k, min(v), k, sd))
head = ("tag", "workload", "ranks", "threads", "ppn")
print("IOK " + " ".join("%s=%s" % (k, best.get(k, "")) for k in head)
      + " reps=%s nok=%d " % (sys.argv[2], len(reps)) + " ".join(stat) + " "
      + " ".join("%s=%s" % (k, v) for k, v in best.items() if k not in head))
PY
}

# igloo_run <tag> <workload> <ranks> <threads> <ppn> [core-mask|NONE] [reps] [check|finite|none]
# One IPIN, IAFF and IREP (or IFAIL) line per repetition, one IOK line per point, one ICHECK line when asked.
# Tags are <workload>-<R>x<T>[-<variant>]. Each rank runs under igloo-stamp.sh (launch_s, life_s, teardown_s).
# Environment: IGLOO_WRAP wraps the binary (numactl --interleave=all); IGLOO_OUTPUT=shm puts OUTPUT/ on /dev/shm;
# IGLOO_INI_SET="key=value ..." edits [IGLOO-General] of the run's input.ini; IGLOO_AA=1 runs the twin snapshot
# (the A/A pair); IGLOO_LAUNCH_RES_OPTIONS is handed to mpirun as RES_OPTIONS (README.md, "Launch cost").
igloo_run () {
  local tag=$1 wl=$2 R=$3 T=$4 ppn=$5 mc reps chk
  if [ $# -ge 6 ]; then
    mc=$6; [ -n "$mc" ] || { echo "IFAIL tag=$tag reason=empty-core-mask"; return 1; }
  else mc=NONE; fi
  reps=${7:-1}; chk=${8:-none}
  local bin=$BIN; [ -n "${IGLOO_AA:-}" ] && bin=$BIN_AA
  [ -f "$CASE/$wl/md5.txt" ] || { echo "IFAIL tag=$tag reason=workload-not-built"; return 1; }
  local wlmd5; wlmd5=$(md5sum < "$CASE/$wl/md5.txt" | cut -c1-12)
  local repfile=$W/work/$SLURM_JOB_ID/$tag.reps
  : > "$repfile"
  local rep D C shm rc t0 t1 wall stamps nstop nrep nwit rsspid mpid opid kept=""
  local -a lenv
  for rep in $(seq 1 "$reps"); do
    D=$W/work/$SLURM_JOB_ID/$tag-r$rep
    igloo_casedir "$D" "$wl" || return 1
    C=$CASEDIR
    shm=""
    if [ "${IGLOO_OUTPUT:-}" = shm ]; then
      shm=/dev/shm/${USER:-igloo}/$SLURM_JOB_ID/$tag-r$rep
      rm -rf "$shm" && mkdir -p "$shm" && rmdir "$C/OUTPUT" && ln -s "$shm" "$C/OUTPUT"
    fi
    [ -n "${IGLOO_INI_SET:-}" ] && igloo_ini_set "$C/input.ini" $IGLOO_INI_SET
    cd "$C" || return 1

    # resident memory every 5 s: the largest rank and the node total of the IGLOO processes
    ( maxrss=0; maxsum=0; me=$(id -u)
      while :; do
        rss=$(ps -o rss=,uid= -C IGLOO 2>/dev/null | awk -v u="$me" '$2==u {print $1}')
        if [ -n "$rss" ]; then
          hi=$(echo "$rss" | sort -n | tail -1)
          sum=$(echo "$rss" | awk '{s+=$1} END{print s+0}')
          [ "$hi" -gt "$maxrss" ] 2>/dev/null && maxrss=$hi
          [ "$sum" -gt "$maxsum" ] 2>/dev/null && maxsum=$sum
        fi
        echo "$maxrss $maxsum" > rss.txt
        sleep 5
      done ) &
    rsspid=$!

    export OMP_NUM_THREADS=$T
    if [ "$mc" = NONE ]; then
      export I_MPI_PIN_DOMAIN=omp
    else
      I_MPI_PIN_DOMAIN=$(igloo_domains "$mc" "$R" "$T") || {
        echo "IFAIL tag=$tag rep=$rep reason=domain-build"; kill $rsspid 2>/dev/null; cd "$W"; continue; }
      export I_MPI_PIN_DOMAIN
    fi
    echo "IPIN tag=$tag rep=$rep I_MPI_PIN_DOMAIN=$I_MPI_PIN_DOMAIN"
    igloo_wait_clear

    lenv=(I_MPI_DEBUG=4 OMP_DISPLAY_AFFINITY=TRUE "IGLOO_STAMPS=$D/stamps")
    [ -n "${IGLOO_LAUNCH_RES_OPTIONS:-}" ] && lenv+=("RES_OPTIONS=$IGLOO_LAUNCH_RES_OPTIONS")
    mkdir -p "$D/stamps"
    t0=$(date +%s%N)
    env "${lenv[@]}" mpirun -n "$R" -ppn "$ppn" bash "$STAMP" ${IGLOO_WRAP:-} "$bin" \
      > run_out.txt 2> run_err.txt < /dev/null &
    mpid=$!
    ( sleep 30; kill -0 $mpid 2>/dev/null && igloo_occupancy "$tag" ) &
    opid=$!
    wait $mpid; rc=$?
    t1=$(date +%s%N)
    kill $opid 2>/dev/null; wait $opid 2>/dev/null
    kill $rsspid 2>/dev/null; wait $rsspid 2>/dev/null
    wall=$(( (t1 - t0) / 1000000 )); wall=$(printf '%d.%03d' $((wall / 1000)) $((wall % 1000)))
    stamps=$(igloo_stamps "$D/stamps" "$t0" "$t1")
    local rss_rank rss_node
    rss_rank=$(awk '{printf "%.3f", $1/1048576}' rss.txt 2>/dev/null); rss_rank=${rss_rank:-0.000}
    rss_node=$(awk '{printf "%.3f", $2/1048576}' rss.txt 2>/dev/null); rss_node=${rss_node:-0.000}
    { cat run_out.txt; echo "---- stderr ----"; cat run_err.txt; } > "$LOGS/$tag-r$rep.log"

    # Scored only after a clean end: exit 0, the normal end line, the report, and the rank witness (a serial
    # binary under mpirun passes every other test while decomposing nothing).
    nstop=$(grep -c 'Stop condition' run_out.txt)
    nrep=$(grep -c '\[IGLOO-timing\] sweep' run_out.txt)
    nwit=1; [ "$R" -gt 1 ] && nwit=$(grep -cE "MPI ranks = $R( |\$)" run_out.txt)
    local parsed aff fields
    parsed=$(igloo_parse "$tag" "$rep" run_out.txt run_err.txt)
    aff=$(echo "$parsed" | sed -n 1p); fields=$(echo "$parsed" | sed -n 2p)
    echo "$aff"
    if [ "$rc" -ne 0 ] || [ "$nstop" -eq 0 ] || [ "$nrep" -eq 0 ] || [ "$nwit" -eq 0 ]; then
      echo "IFAIL tag=$tag rep=$rep rc=$rc stop=$nstop report=$nrep rank_witness=$nwit wall_s=$wall $stamps" \
           "-- log $LOGS/$tag-r$rep.log"
      tail -12 run_out.txt | sed 's/^/  | /'; tail -6 run_err.txt | sed 's/^/  ! /'
      cd "$W"; rm -rf "$D"; igloo_rmshm "$shm"
      continue
    fi
    echo "IREP tag=$tag workload=$wl ranks=$R threads=$T ppn=$ppn rep=$rep rc=$rc wall_s=$wall $stamps" \
         "rss_rank_gb=$rss_rank rss_node_gb=$rss_node wl_md5=$wlmd5 $fields" | tee -a "$repfile"
    cd "$W"
    if [ "$rep" -lt "$reps" ] || [ "$chk" = none ]; then
      rm -rf "$D"; igloo_rmshm "$shm"
    else
      kept=$D
    fi
  done
  igloo_summary "$repfile" "$reps"

  # the workload's own oracles on the last repetition's output
  if [ "$chk" != none ]; then
    if [ -z "$kept" ]; then
      echo "ICHECK tag=$tag finite=skip check=skip reason=last-repetition-failed"
    else
      local frc crc=skip
      cd "$C"
      python3 -B "$SRC/tests/tools/check_finite.py" > "$LOGS/$tag-finite.log" 2>&1; frc=$?
      if [ "$chk" = check ] && [ -f check.py ]; then
        IGLOO_EXE=$bin python3 -B check.py > "$LOGS/$tag-check.log" 2>&1; crc=$?
      fi
      echo "ICHECK tag=$tag finite=$frc check=$crc"
      tail -4 "$LOGS/$tag-finite.log" | sed 's/^/  | /'
      [ -f "$LOGS/$tag-check.log" ] && tail -4 "$LOGS/$tag-check.log" | sed 's/^/  | /'
      cd "$W"; rm -rf "$kept"; igloo_rmshm "$shm"
    fi
  fi
  return 0
}

# igloo_matrix <workload> <reps> <check|finite>: one workload's scaling curve on the job's node -- the 1x1 anchor
# (with the workload's checks), every split of 4, 8, 16 and 32 slots spread over the sockets, the 80-slot splits
# on Intel MPI's omp domains (the box without 80 x 1: 80 ranks at the finalize peak exceed the node's memory), and
# the A/A twin at 1x1 (not the box) and 4x20.
igloo_matrix () {
  local wl=$1 reps=$2 chk=$3 N R Tn RT
  echo "=== $wl: the 1x1 anchor ==="
  igloo_run "$wl-1x1" "$wl" 1 1 1 "$(igloo_spread 1)" "$reps" "$chk"
  echo "=== $wl: every split of 4, 8, 16 and 32 slots, spread over the sockets ==="
  for N in 4 8 16 32; do
    for R in 1 2 4 8 16 32; do
      [ "$R" -le "$N" ] || continue
      Tn=$(( N / R ))
      igloo_run "$wl-${R}x${Tn}" "$wl" "$R" "$Tn" "$R" "$(igloo_spread "$N")" "$reps"
    done
  done
  echo "=== $wl: the 80-slot splits ==="
  for RT in 1x80 4x20 8x10 16x5 20x4 40x2 80x1; do
    [ "$wl" = box-4fam ] && [ "$RT" = 80x1 ] && continue
    R=${RT%x*}; Tn=${RT#*x}
    igloo_run "$wl-$RT" "$wl" "$R" "$Tn" "$R" NONE "$reps"
  done
  echo "=== $wl: the A/A twin, a second snapshot of the same binary ==="
  [ "$wl" = box-4fam ] || IGLOO_AA=1 igloo_run "$wl-1x1-aa" "$wl" 1 1 1 "$(igloo_spread 1)" "$reps"
  IGLOO_AA=1 igloo_run "$wl-4x20-aa" "$wl" 4 20 4 NONE "$reps"
}
