#!/bin/bash
#SBATCH --job-name=igloo-io
#SBATCH --partition=free
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=80
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=03:00:00
#
# Where the box's gas read and grid write spend their time, on one exclusive node. IGLOO_IO_STEPS picks the steps
# (blank-, comma- or colon-separated; default all, in this order):
#   facts     the file systems of the harness root, the job's directory and /dev/shm, the dirty-page limits, the
#             tools, the compiler and the default of its buffered_io option
#   write     the solver's grid write: box-4fam at 1 x 80 as is, under FORT_BUFFERED=true, with OUTPUT/ on tmpfs,
#             both (IGLOO_REPS repetitions, default 3); khrt-stress at 1 x 32 and 8 x 10 as is and buffered
#             (IGLOO_IO_KHRT_REPS, default 5); every run under io-times.sh (its user and system time: ITIMES lines)
#   read      the probe on the box at one thread: ORION's reader, one formatted line read per line, stream reads
#             (cold and warm), the tokenizer and six conversion kernels; the random sets of make_values.py; the
#             suite's .tec files
#   threads   the probe's kernels at 1, 2, 4, 10 and 20 threads on one socket and at 40 and 80 spread
#   stream    a 3.19 GB buffer written in 64 MiB blocks to the job's directory and to /dev/shm
#   bcast     the gas broadcast of the box's field at 4, 8, 16 and 40 ranks of one thread
#   buffered  ORION's writer on 10 M values with -assume buffered_io on the main program, on the module holding
#             the OPEN, on both, and under FORT_BUFFERED=true; the write calls counted by strace
# The solver points run $IGLOO_PERF_BIN (default $IGLOO_PERF_ROOT/bin/IGLOO-hyb); the probe's executables come
# from $IGLOO_IO_PROBE_DIR (default $IGLOO_PERF_ROOT/bin), installed by `build-cfg.sh probe`. Submit from the
# harness root:
#   sbatch --output=$IGLOO_PERF_ROOT/logs/%x-%j.out --export=ALL,IGLOO_PERF_ROOT=$IGLOO_PERF_ROOT scripts/igloo-io.slurm
set -u
W=${IGLOO_PERF_ROOT:?}/io
source "$IGLOO_PERF_ROOT/scripts/_igloo_common.sh"
igloo_setup || exit 1
STEPS=$(echo "${IGLOO_IO_STEPS:-facts write read threads stream bcast buffered}" | tr ',:' '  ')
REPS=${IGLOO_REPS:-3}
KREPS=${IGLOO_IO_KHRT_REPS:-5}
TW="bash $ROOT/scripts/io-times.sh"
JW=$W/work/$SLURM_JOB_ID
SHM=/dev/shm/${USER:-igloo}/$SLURM_JOB_ID
BOX=$CASE/box-4fam/INPUT/solfile.tec
echo "STEPS $STEPS"

# The ITIMES lines of a point's repetitions, from their kept logs, tagged with the point and the repetition.
io_times () {
  local f r
  for f in "$LOGS/$1"-r*.log; do
    [ -f "$f" ] || continue
    r=${f##*-r}; r=${r%.log}
    grep '^ITIMES ' "$f" | sed "s/^ITIMES /ITIMES tag=$1 rep=$r /"
  done
}

io_facts () {
  local d man
  for d in "$ROOT" "$JW" /dev/shm; do
    mkdir -p "$d"
    echo "IFACT fs path=$d type=$(stat -f -c %T "$d") $(df -PT "$d" | awk 'NR == 2 {print "source=" $1, "fstype=" $2, "size_kb=" $3, "avail_kb=" $5}')"
    findmnt -n -o TARGET,FSTYPE,OPTIONS -T "$d" 2>/dev/null | sed 's/^/IFACT mount /'
  done
  echo "IFACT vm dirty_ratio=$(cat /proc/sys/vm/dirty_ratio) dirty_background_ratio=$(cat /proc/sys/vm/dirty_background_ratio)" \
       "dirty_bytes=$(cat /proc/sys/vm/dirty_bytes) dirty_background_bytes=$(cat /proc/sys/vm/dirty_background_bytes)" \
       "dirty_expire_centisecs=$(cat /proc/sys/vm/dirty_expire_centisecs)"
  echo "IFACT tools strace=$(command -v strace || echo none) time=$([ -x /usr/bin/time ] && echo /usr/bin/time || echo none)" \
       "numactl=$(command -v numactl || echo none) dd=$(dd --version 2>/dev/null | head -1 | tr ' ' '_') kernel=$(uname -r)"
  echo "IFACT compiler $(ifx --version 2>&1 | head -1)"
  man=$(ls "$(dirname "$(command -v ifx)")"/../../documentation/en/man/common/man1/ifx.1 2>/dev/null | head -1)
  if [ -n "$man" ]; then
    awk '/^nobuffered_io/ {n = 12} n > 0 {print "IFACT ifx.1 " $0; n--}' "$man"
  else
    echo "IFACT ifx.1 not found"
  fi
  ldd "$BIN" | grep -E 'ifcore|ifport|imf|svml|iomp|libmpi' | sed 's/^[[:space:]]*/IFACT ldd /'
  lscpu -p=CPU,CORE,SOCKET,NODE | grep -v '^#' | tr '\n' ' ' | sed 's/^/IFACT cpus /; s/$/\n/'
}

io_write () {
  echo "=== the grid write: box-4fam at 1 x 80, as is / buffered / OUTPUT on tmpfs / both ==="
  IGLOO_WRAP="$TW" igloo_run box-4fam-1x80-unbuf box-4fam 1 80 1 NONE "$REPS"
  io_times box-4fam-1x80-unbuf
  IGLOO_WRAP="$TW env FORT_BUFFERED=true" igloo_run box-4fam-1x80-buf box-4fam 1 80 1 NONE "$REPS" finite
  io_times box-4fam-1x80-buf
  IGLOO_OUTPUT=shm IGLOO_WRAP="$TW" igloo_run box-4fam-1x80-unbuf-shm box-4fam 1 80 1 NONE "$REPS"
  io_times box-4fam-1x80-unbuf-shm
  IGLOO_OUTPUT=shm IGLOO_WRAP="$TW env FORT_BUFFERED=true" igloo_run box-4fam-1x80-buf-shm box-4fam 1 80 1 NONE "$REPS"
  io_times box-4fam-1x80-buf-shm
  echo "=== the particle streams: khrt-stress at 1 x 32 and 8 x 10, as is / buffered ==="
  IGLOO_WRAP="$TW" igloo_run khrt-stress-1x32-unbuf khrt-stress 1 32 1 "$(igloo_spread 32)" "$KREPS" check
  io_times khrt-stress-1x32-unbuf
  IGLOO_WRAP="$TW env FORT_BUFFERED=true" igloo_run khrt-stress-1x32-buf khrt-stress 1 32 1 "$(igloo_spread 32)" "$KREPS" check
  io_times khrt-stress-1x32-buf
  IGLOO_WRAP="$TW" igloo_run khrt-stress-8x10-unbuf khrt-stress 8 10 8 NONE "$KREPS" check
  io_times khrt-stress-8x10-unbuf
  IGLOO_WRAP="$TW env FORT_BUFFERED=true" igloo_run khrt-stress-8x10-buf khrt-stress 8 10 8 NONE "$KREPS" check
  io_times khrt-stress-8x10-buf
}

# The probe's executables, copied into the job's directory (a rebuild during the job cannot reach them).
PROBE=
io_probe_setup () {
  local d=${IGLOO_IO_PROBE_DIR:-$ROOT/bin} e
  [ -n "$PROBE" ] && return 0
  mkdir -p "$W/bin.$SLURM_JOB_ID/probe"
  for e in io-probe io-probe-wr io-probe-wr-main io-probe-wr-mod io-probe-wr-both; do
    [ -x "$d/$e" ] || { echo "IFAIL probe reason=no-executable ($d/$e)"; return 1; }
    cp -f "$d/$e" "$W/bin.$SLURM_JOB_ID/probe/$e"
    echo "PROBE $e md5=$(md5sum < "$W/bin.$SLURM_JOB_ID/probe/$e" | cut -d' ' -f1)"
  done
  if [ -f "$d/io-probe.commit" ]; then echo "PROBE-COMMIT $(tr '\n' ' ' < "$d/io-probe.commit")"
  else echo "PROBE-COMMIT unknown (no $d/io-probe.commit)"; fi
  PROBE=$W/bin.$SLURM_JOB_ID/probe
}

io_read () {
  local list=$JW/suite-tec.txt
  echo "=== the read split on the box, one thread ==="
  dd if="$BOX" of=/dev/null bs=64M iflag=direct 2>&1 | tail -1 | sed 's/^/IFACT dd-direct /'
  OMP_NUM_THREADS=1 "$PROBE/io-probe" read-split "$BOX"
  echo "=== the random sets, one thread ==="
  mkdir -p "$SHM/values"
  python3 -B "$SRC/tests/bench/monolith/io-probe/make_values.py" "$SHM/values" 10000000
  OMP_NUM_THREADS=1 "$PROBE/io-probe" randsets "$SHM/values"
  rm -rf "$SHM/values"
  echo "=== the suite's .tec files, one thread ==="
  find "$SRC/tests" -name '*.tec' ! -path '*/OUTPUT/*' ! -path '*/ref-*' -exec realpath {} \; | LC_ALL=C sort -u > "$list"
  echo "IFACT suite-tec files=$(wc -l < "$list") bytes=$(xargs cat < "$list" | wc -c)"
  OMP_NUM_THREADS=1 "$PROBE/io-probe" census "$list"
}

# Every core once, socket after socket (the first CPU of each core): the places of the thread-scaling step.
io_places () {
  python3 - <<'PY'
import subprocess
out = subprocess.run(["lscpu", "-p=CPU,CORE,SOCKET"], stdout=subprocess.PIPE).stdout.decode()
first = {}
for line in out.splitlines():
    if line.startswith("#") or not line.strip():
        continue
    cpu, core, s = (int(x) for x in line.split(",")[:3])
    if (s, core) not in first or cpu < first[(s, core)]:
        first[(s, core)] = cpu
print(",".join("{%d}" % first[k] for k in sorted(first)))
PY
}

io_threads () {
  local places; places=$(io_places)
  echo "=== the kernels at 1-20 threads on one socket, 40 and 80 spread ==="
  echo "IFACT places $places"
  OMP_PLACES="$places" OMP_PROC_BIND=close OMP_NUM_THREADS=80 \
    "$PROBE/io-probe" scale "$BOX" 0 1c 2c 4c 10c 20c 40s 80s
}

io_stream () {
  local f rep
  echo "=== the stream write: 3.19 GB in 64 MiB blocks ==="
  mkdir -p "$SHM"
  for f in "$JW/stream.bin" "$SHM/stream.bin"; do
    for rep in 1 2; do
      "$PROBE/io-probe" swrite "$f" 3193380360 67108864
      rm -f "$f"
    done
  done
}

io_bcast () {
  local R mask
  local -a lenv
  echo "=== the gas broadcast: the box's field at 4, 8, 16 and 40 ranks of one thread ==="
  for R in 4 8 16 40; do
    mask=$(igloo_domains "$(igloo_spread "$R")" "$R" 1) || { echo "IFAIL bcast ranks=$R reason=domain-build"; continue; }
    lenv=(OMP_NUM_THREADS=1 "I_MPI_PIN_DOMAIN=$mask" I_MPI_DEBUG=4)
    [ -n "${IGLOO_LAUNCH_RES_OPTIONS:-}" ] && lenv+=("RES_OPTIONS=$IGLOO_LAUNCH_RES_OPTIONS")
    env "${lenv[@]}" mpirun -n "$R" -ppn "$R" "$PROBE/io-probe" bcast 301 101 101 18 5 < /dev/null
    echo "IBCAST ranks=$R rc=$?"
  done
}

# io_strace <label> <command ...>: the write system calls of one writer run (strace -c), parsed and raw.
io_strace () {
  local label=$1 s=$JW/strace-$1.txt
  shift
  strace -f -c -e trace=write,pwrite64,writev -o "$s" "$@" > /dev/null
  awk -v l="$label" '$NF ~ /^(write|pwrite64|writev)$/ {print "ISTRACE variant=" l " syscall=" $NF " calls=" $4 " seconds=" $2}' "$s"
  sed 's/^/ISTRACE-RAW /' "$s"
  rm -f "$JW/wr.tec"
}

io_buffered () {
  local f=$JW/wr.tec e
  echo "=== -assume buffered_io: main program / module holding the OPEN / both; FORT_BUFFERED ==="
  for e in io-probe-wr io-probe-wr-main io-probe-wr-mod io-probe-wr-both; do
    $TW "$PROBE/$e" "$f" 100 2>&1; rm -f "$f"
  done
  FORT_BUFFERED=true $TW "$PROBE/io-probe-wr" "$f" 100 2>&1; rm -f "$f"
  if command -v strace > /dev/null; then
    for e in io-probe-wr io-probe-wr-main io-probe-wr-mod io-probe-wr-both; do
      io_strace "$e" "$PROBE/$e" "$f" 48
    done
    io_strace io-probe-wr-FORT_BUFFERED env FORT_BUFFERED=true "$PROBE/io-probe-wr" "$f" 48
  fi
}

WL=""
for s in $STEPS; do
  case $s in
    write) WL="$WL box-4fam khrt-stress" ;;
    read|threads) WL="$WL box-4fam" ;;
  esac
done
WL=$(echo $WL | tr ' ' '\n' | awk 'NF && !seen[$0]++' | tr '\n' ' ')
if [ -n "$WL" ]; then
  bash "$ROOT/scripts/make-workloads.sh" $WL || { echo "IGLOO-IO STOPPED: the workloads did not build"; exit 1; }
fi
mkdir -p "$JW"

for s in $STEPS; do
  echo "##### step $s $(date -u +%Y-%m-%dT%H:%M:%SZ)"
  case $s in
    facts) io_facts ;;
    write) io_write ;;
    read)     io_probe_setup && io_read ;;
    threads)  io_probe_setup && io_threads ;;
    stream)   io_probe_setup && io_stream ;;
    bcast)    io_probe_setup && io_bcast ;;
    buffered) io_probe_setup && io_buffered ;;
    *) echo "IFAIL step=$s reason=unknown-step" ;;
  esac
done
rm -rf "$SHM"
echo "IGLOO-IO DONE $(date -u +%Y-%m-%dT%H:%M:%SZ)"
