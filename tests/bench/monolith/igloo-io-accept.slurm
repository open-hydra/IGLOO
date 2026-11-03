#!/bin/bash
#SBATCH --job-name=igloo-io-accept
#SBATCH --partition=free
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=80
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=04:00:00
#
# The box's gas read and grid write at the 1x1 anchor (with check_finite.py) and at 1x80, 4x20, 40x2 and the A/A
# twin at 4x20 (IGLOO_ACCEPT_POINTS: <R>x<T>, -aa the twin, -shm OUTPUT/ on /dev/shm), IGLOO_REPS repetitions
# (default 3), with the file system's floors measured right before every repetition (one IFLOOR line):
#   F  a buffer of the box's grid output size (3.19 GB) written in 64 MiB blocks into the repetition's case
#      directory, timed to the end of close, then deleted (io-probe swrite)
#   T  the same into /dev/shm
#   C  the box's gas file read whole into a fresh buffer in 256 MiB reads, timed from open to close, with the
#      node's page cache dropped and the file server's warm: its pages dropped (posix_fadvise), the file read once
#      untimed (C_warm_s), dropped, the timed read, dropped again; fincore counts the resident bytes after the last
#      three steps (evict_ok=1: at most a thousandth of the file resident after each drop)
# Then the outputs of one split with one thread per rank (byte-reproducible: the cross-rank sum runs in rank order)
# from a second binary and from the job's, compared file by file: one ICOMPARE line, the two runs' lines in
# logs/<job>/compare.out. IGLOO_ACCEPT_COMPARE="<workload> <ranks>" (default "box-4fam 4"; "none" skips it),
# IGLOO_ACCEPT_BASE the second binary (default $IGLOO_PERF_ROOT/bin-base/IGLOO-hyb: `build-cfg.sh hyb <tree>
# $IGLOO_PERF_ROOT/bin-base`); the probe from $IGLOO_IO_PROBE_DIR (default $IGLOO_PERF_ROOT/bin, `build-cfg.sh
# probe`). Submit from the harness root:
#   sbatch --output=$IGLOO_PERF_ROOT/logs/%x-%j.out --export=ALL,IGLOO_PERF_ROOT=$IGLOO_PERF_ROOT \
#     scripts/igloo-io-accept.slurm
set -u
W=${IGLOO_PERF_ROOT:?}/accept
source "$IGLOO_PERF_ROOT/scripts/_igloo_common.sh"
igloo_setup || exit 1
echo "HARNESS $(basename "$0") md5=$(md5sum < "$0" | cut -d' ' -f1) _igloo_common.sh" \
     "md5=$(md5sum < "$ROOT/scripts/_igloo_common.sh" | cut -d' ' -f1)"
POINTS=${IGLOO_ACCEPT_POINTS:-1x1 1x80 4x20 40x2 4x20-aa}
REPS=${IGLOO_REPS:-3}
read -r CWL CR <<< "${IGLOO_ACCEPT_COMPARE:-box-4fam 4}"
BOX=$CASE/box-4fam/INPUT/solfile.tec
FBYTES=3193380360
FBLOCK=67108864
TDIR=/dev/shm/${USER:-igloo}/$SLURM_JOB_ID

WL=box-4fam
[ "$CWL" = none ] || [ "$CWL" = box-4fam ] || WL="$WL $CWL"
bash "$ROOT/scripts/make-workloads.sh" $WL || { echo "IGLOO-ACCEPT STOPPED: the workloads did not build"; exit 1; }

# The probe and the second binary, copied into the job's directory (a rebuild during the job cannot reach them).
PD=${IGLOO_IO_PROBE_DIR:-$ROOT/bin}
[ -x "$PD/io-probe" ] || { echo "IGLOO-ACCEPT STOPPED: no $PD/io-probe"; exit 1; }
mkdir -p "$W/bin.$SLURM_JOB_ID/probe"
PROBE=$W/bin.$SLURM_JOB_ID/probe/io-probe
cp -f "$PD/io-probe" "$PROBE"
echo "PROBE io-probe md5=$(md5sum < "$PROBE" | cut -d' ' -f1)"
if [ -f "$PD/io-probe.commit" ]; then echo "PROBE-COMMIT $(tr '\n' ' ' < "$PD/io-probe.commit")"
else echo "PROBE-COMMIT unknown (no $PD/io-probe.commit)"; fi
BASE=
if [ "$CWL" != none ]; then
  B0=${IGLOO_ACCEPT_BASE:-$ROOT/bin-base/IGLOO-hyb}
  [ -x "$B0" ] || { echo "IGLOO-ACCEPT STOPPED: no second binary $B0"; exit 1; }
  mkdir -p "$W/bin.$SLURM_JOB_ID-base"
  BASE=$W/bin.$SLURM_JOB_ID-base/IGLOO
  cp -f "$B0" "$BASE"
  echo "BASE-BINARY md5=$(md5sum < "$BASE" | cut -d' ' -f1) $(stat -c %s "$BASE") bytes (snapshot: $BASE)"
  if [ -f "$B0.commit" ]; then echo "BASE-COMMIT $(tr '\n' ' ' < "$B0.commit")"
  else echo "BASE-COMMIT unknown (no $B0.commit)"; fi
fi
echo "FLOORS F and T bytes=$FBYTES block=$FBLOCK into each repetition's case directory and $TDIR;" \
     "C file=$BOX bytes=$(stat -c %s "$BOX")"

# accept_floors <tag> <rep> <case directory>: the IFLOOR line of one repetition (IGLOO_PRE_REP).
accept_floors () {
  mkdir -p "$TDIR"
  python3 - "$1" "$2" "$PROBE" "$3/floor-F.bin" "$TDIR/floor-T.bin" "$FBYTES" "$FBLOCK" "$BOX" <<'PY'
import mmap, os, re, subprocess, sys, time
tag, rep, probe, ffile, tfile, fbytes, fblock, box = sys.argv[1:9]


def swrite(path):
    p = subprocess.run([probe, "swrite", path, fbytes, fblock], stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    line = next((l for l in p.stdout.decode(errors="replace").splitlines() if l.startswith("IPROBE step=stream")), "")
    try:
        os.remove(path)
    except OSError:
        pass
    return dict(re.findall(r"(\w+)=(\S+)", line)), p.returncode


def sec(kv, k):
    try:
        return "%.3f" % float(kv[k])
    except (KeyError, ValueError):
        return "nan"


def resident():
    r = subprocess.run(["fincore", "-b", "-n", "-o", "RES", box], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    try:
        return int(r.stdout.split()[0])
    except (IndexError, ValueError):
        return -1


def drop():
    fd = os.open(box, os.O_RDONLY)
    try:
        os.posix_fadvise(fd, 0, 0, os.POSIX_FADV_DONTNEED)
    finally:
        os.close(fd)


def read_whole(fresh):
    """The file in 256 MiB reads, open to close: into a fresh file-sized buffer (its pages faulted in by the reads),
    or through one 256 MiB buffer."""
    t0 = time.perf_counter()
    fd = os.open(box, os.O_RDONLY)
    n = os.fstat(fd).st_size
    buf = mmap.mmap(-1, n if fresh else min(n, 256 << 20))
    mv = memoryview(buf)
    got = 0
    while got < n:
        a = got if fresh else 0
        k = os.readv(fd, [mv[a:a + min(256 << 20, n - got)]])
        if k <= 0:
            break
        got += k
    os.close(fd)
    t1 = time.perf_counter()
    mv.release()
    buf.close()
    return t1 - t0, got, n


fkv, frc = swrite(ffile)
tkv, trc = swrite(tfile)
drop()
tw, gw, n = read_whole(False)
drop()
r0 = resident()
tc, got, n = read_whole(True)
r1 = resident()
drop()
r2 = resident()
lim = n // 1000
ok = int(0 <= r0 <= lim and 0 <= r2 <= lim and got == n and gw == n)
print("IFLOOR tag=%s rep=%s F_s=%s F_write_s=%s F_close_s=%s F_rc=%d T_s=%s T_rc=%d C_warm_s=%.3f C_s=%.3f C_bytes=%d "
      "C_MBps=%.0f res_dropped=%d res_read=%d res_after=%d evict_ok=%d" % (
          tag, rep, sec(fkv, "wall_s"), sec(fkv, "write_s"), sec(fkv, "close_s"), frc, sec(tkv, "wall_s"), trc, tw,
          tc, got, got / 1e6 / tc, r0, r1, r2, ok))
PY
}

# accept_md5 <tag> <rep> <case directory>: the md5 of every output file (IGLOO_POST_REP of the comparison).
accept_md5 () {
  ( cd "$3/OUTPUT" && find . -type f | LC_ALL=C sort | xargs -r md5sum ) > "$LOGS/$1-r$2.md5"
}

NC=$(lscpu -p=CORE,SOCKET | grep -v '^#' | sort -u | wc -l)
for pt in $POINTS; do
  [[ $pt =~ ^[0-9]+x[0-9]+(-aa|-shm)?$ ]] || { echo "IFAIL point=$pt reason=not-RxT[-aa|-shm]"; continue; }
  RT=${pt%%-*}; R=${RT%x*}; T=${RT#*x}
  mc=NONE; [ $((R * T)) -lt "$NC" ] && mc=$(igloo_spread $((R * T)))
  chk=none; [ "$pt" = 1x1 ] && chk=finite
  echo "=== box-4fam-$pt ==="
  case $pt in
    *-aa) IGLOO_PRE_REP=accept_floors IGLOO_AA=1 igloo_run "box-4fam-$pt" box-4fam "$R" "$T" "$R" "$mc" "$REPS" "$chk" ;;
    *-shm) IGLOO_PRE_REP=accept_floors IGLOO_OUTPUT=shm igloo_run "box-4fam-$pt" box-4fam "$R" "$T" "$R" "$mc" "$REPS" "$chk" ;;
    *) IGLOO_PRE_REP=accept_floors igloo_run "box-4fam-$pt" box-4fam "$R" "$T" "$R" "$mc" "$REPS" "$chk" ;;
  esac
done
rmdir "$TDIR" 2>/dev/null

if [ "$CWL" != none ]; then
  echo "=== $CWL at ${CR}x1: the second binary, then the job's; every output file compared ==="
  {
    IGLOO_POST_REP=accept_md5 BIN=$BASE igloo_run "$CWL-${CR}x1-base" "$CWL" "$CR" 1 "$CR" "$(igloo_spread "$CR")" 1
    IGLOO_POST_REP=accept_md5 igloo_run "$CWL-${CR}x1" "$CWL" "$CR" 1 "$CR" "$(igloo_spread "$CR")" 1
  } > "$LOGS/compare.out" 2>&1
  python3 - "$LOGS" "$CWL-${CR}x1" "$LOGS/compare.out" "$(md5sum < "$BASE" | cut -d' ' -f1)" <<'PY'
import os, sys
logs, tag, out, bmd5 = sys.argv[1:5]
lines = open(out, errors="replace").read().splitlines()


def ok(t):
    return int(any(l.startswith("IREP tag=%s " % t) for l in lines))


def load(t):
    try:
        return dict(reversed(l.rstrip("\n").split(None, 1)) for l in open(os.path.join(logs, t + "-r1.md5")))
    except (OSError, ValueError):
        return {}


a, b = load(tag + "-base"), load(tag)
common = sorted(set(a) & set(b))
diff = [f for f in common if a[f] != b[f]]
only_a, only_b = sorted(set(a) - set(b)), sorted(set(b) - set(a))
good = ok(tag + "-base") and ok(tag) and a and not diff and not only_a and not only_b
print("ICOMPARE tag=%s base_md5=%s base_ok=%d job_ok=%d files=%d identical=%d differ=%d only_base=%d only_job=%d "
      "verdict=%s" % (tag, bmd5, ok(tag + "-base"), ok(tag), len(common), len(common) - len(diff), len(diff),
                      len(only_a), len(only_b), "identical" if good else "DIFFERENT"))
for f in diff + only_a + only_b:
    print("ICOMPARE-FILE tag=%s file=%s base=%s job=%s" % (tag, f, a.get(f, "-"), b.get(f, "-")))
PY
fi
echo "IGLOO-ACCEPT DONE $(date -u +%Y-%m-%dT%H:%M:%SZ)"
