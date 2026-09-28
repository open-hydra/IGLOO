#!/usr/bin/env python3
"""B-VAL-6 thread-invariance gate: KHRT output must not depend on the thread count.

`obj_IGLOO%solve` runs `integrate` inside `!$OMP PARALLEL DO SCHEDULE(DYNAMIC)`, and the KH
shed path adds per-parent state that is written from inside that region and drained after it.
Two distinct things can go wrong there, and only one of them is visible to the other gates:

  * ORDERING — which child lands on which particle slot. The drain walks `ip` ascending and
    then push order within each list, so the mapping is fixed by the traversal rather than by
    the schedule. A sorted-multiset compare deliberately ABSORBS ordering differences: record
    order in `trajectories-*.dat` is OMP-nondeterministic by design and is not a defect.
  * TEARING — `unitTraj`/`unitExit` are written from inside the parallel region with no
    synchronisation, so two threads can interleave mid-record. A sorted multiset does NOT
    absorb that: a torn line is a different line. This gate is the only thing in the suite
    that would catch it.

Established as a PASSING BASELINE while the one-shed cap was still enforced, so that when the
cap is lifted a red result here means the lift broke thread-invariance, not that it was already
broken. Compares OMP_NUM_THREADS = 1, 2, 4.

Usage: check_threads.py <path-to-IGLOO-binary>
"""
import glob
import hashlib
import os
import re
import resource
import shutil
import subprocess
import sys
import time

THREADS = (1, 2, 4)
FILES = ("OUTPUT/trajectories-A.dat", "OUTPUT/outloc-A.dat")
EXITFILE = "OUTPUT/outloc-A.dat"
ZONE_HEADERS = 2          #> one `variables=` line + one per-group `Zone` line
DIAG = "diag_threads.txt"
SENTINEL = "OUTPUT/.sentinel"


def digest(path):
    """Hash of the record multiset: sorted lines, so ordering is ignored but content is not."""
    with open(path) as fh:
        lines = sorted(fh.readlines())
    h = hashlib.md5()
    for ln in lines:
        h.update(ln.encode())
    return h.hexdigest(), len(lines)


def exe_stat(exe):
    """Identity of the binary at launch time. bin/IGLOO is ONE link target shared by every build
    tree and every e2e case, so a mode/inode/size change between launches is the first thing to
    rule out when this gate misbehaves (BUGS.md O13)."""
    try:
        st = os.stat(exe)
        return (f"mode={oct(st.st_mode)} size={st.st_size} ino={st.st_ino} "
                f"nlink={st.st_nlink} mtime={st.st_mtime_ns}")
    except OSError as exc:
        return f"stat failed: {exc}"


def diag(msg):
    """Failure-path only: never called on a passing run, so PASS stdout stays byte-stable."""
    print(f"  [diag] {msg}")
    with open(DIAG, "a") as fh:
        fh.write(f"{time.time():.3f} {msg}\n")


def dump_context(exe, nthreads, rc):
    diag(f"context at OMP_NUM_THREADS={nthreads}, rc={rc}")
    diag(f"exe now: {exe_stat(exe)}")
    diag(f"cwd: {os.getcwd()}")
    for p in sorted(glob.glob("OUTPUT/*")):
        try:
            diag(f"  OUTPUT: {p} bytes={os.path.getsize(p)}")
        except OSError as exc:
            diag(f"  OUTPUT: {p} unreadable: {exc}")
    for p in ("run_out.txt", "run_err.txt"):
        try:
            tail = open(p, errors="replace").read()[-400:]
            diag(f"  {p} ({os.path.getsize(p)} bytes) tail: {tail!r}")
        except OSError as exc:
            diag(f"  {p} unreadable: {exc}")


def unlimited_stack():
    """Every other e2e gate runs the solver under `bash -c "ulimit -s unlimited"`; this one is
    invoked as a bare python3 command and so inherited the default 8 MB. The project's own runtime
    wrapper treats unlimited as required, so match it rather than run the solver under conditions
    no other gate uses. Best-effort: a refused raise is not a reason to fail the gate."""
    try:
        soft, hard = resource.getrlimit(resource.RLIMIT_STACK)
        if soft != resource.RLIM_INFINITY:
            resource.setrlimit(resource.RLIMIT_STACK, (hard, hard))
    except (ValueError, OSError):
        pass


def run(exe, nthreads):
    shutil.rmtree("OUTPUT", ignore_errors=True)
    os.makedirs("OUTPUT", exist_ok=True)
    open(SENTINEL, "w").close()
    env = dict(os.environ, OMP_NUM_THREADS=str(nthreads), KMP_STACKSIZE="100M")
    before = exe_stat(exe)
    #> Trap the launch itself: an OSError here (O13 shape A was EACCES on execve of the shared
    #  bin/IGLOO) used to surface as a bare traceback with no evidence attached.
    try:
        with open("run_out.txt", "w") as out, open("run_err.txt", "w") as err:
            rc = subprocess.call([exe], stdout=out, stderr=err, env=env)
    except OSError as exc:
        print(f"[FAIL] could not launch solver at OMP_NUM_THREADS={nthreads}: {exc}")
        diag(f"launch OSError errno={exc.errno} ({exc.strerror}) exe={exc.filename}")
        diag(f"exe before launch: {before}")
        dump_context(exe, nthreads, "launch-failed")
        return None
    if rc >= 128:
        print(f"[FAIL] solver died on signal (exit {rc}) at OMP_NUM_THREADS={nthreads}")
        diag(f"exe before launch: {before}")
        dump_context(exe, nthreads, rc)
        return None
    if rc != 0:
        print(f"[FAIL] solver exit {rc} at OMP_NUM_THREADS={nthreads}")
        diag(f"exe before launch: {before}")
        dump_context(exe, nthreads, rc)
        return None
    log = open("run_out.txt", errors="replace").read()
    nkids = sum(int(m) for m in re.findall(r"number of children =\s*(\d+)", log))
    nparents = sum(int(m) for m in re.findall(r"number of particles =\s*(\d+)", log))
    try:
        res = {f: digest(f) for f in FILES}
    except FileNotFoundError as exc:
        print(f"[FAIL] {exc.filename} missing at OMP_NUM_THREADS={nthreads}")
        diag(f"exe before launch: {before}")
        dump_context(exe, nthreads, rc)
        return None

    #> ABSOLUTE completeness check -- the rest of this gate is purely RELATIVE (1 vs 2 vs 4
    #  threads), so a truncation that hit all three runs equally would pass it silently. Every
    #  parcel leaves the domain exactly once and writes exactly one unitExit record, so
    #  exits == parents + children is an exact, build-independent invariant (verified 266 = 25+241
    #  here and 726 = 25+701 on khrt-stress). This is what makes BUGS.md O13 shape B -- exit 0,
    #  empty stderr, a fraction of the records -- fail as TRUNCATION instead of masquerading as a
    #  thread-invariance violation, or passing.
    nexit = res[EXITFILE][1] - ZONE_HEADERS
    if nparents > 0 and nexit != nparents + nkids:
        print(f"[FAIL] incomplete output at OMP_NUM_THREADS={nthreads}: {nexit} exit records "
              f"but {nparents} parents + {nkids} children = {nparents + nkids} parcels")
        print("  [diag] every parcel exits exactly once, so this run did NOT produce complete "
              "output -- suspect BUGS.md O13, not the shed drain.")
        diag(f"exe before launch: {before}")
        diag(f"exits={nexit} parents={nparents} children={nkids}")
        dump_context(exe, nthreads, rc)
        return None

    #> Sentinel: proves whether a FOREIGN process wiped this OUTPUT/ mid-run. The three tests
    #  sharing this directory are serialised by ctest's RESOURCE_LOCK, so a missing sentinel would
    #  mean that assumption is false -- the single most useful fact for closing O13.
    if not os.path.exists(SENTINEL):
        print(f"[FAIL] OUTPUT/ was wiped by another process during the run "
              f"(OMP_NUM_THREADS={nthreads})")
        diag("SENTINEL GONE -- a foreign process removed OUTPUT/ mid-run; the RESOURCE_LOCK "
             "assumption in BUGS.md O13 is refuted")
        dump_context(exe, nthreads, rc)
        return None

    return res, nkids, before


def main():
    if len(sys.argv) < 2:
        print("[FAIL] usage: check_threads.py <path-to-IGLOO-binary>")
        return 1
    exe = sys.argv[1]

    print("KHRT thread-invariance gate (sorted-multiset compare):")
    print(f"{'threads':>7} {'children':>8}  " +
          "  ".join(f"{os.path.basename(f):>24}" for f in FILES))

    if os.path.exists(DIAG):
        os.remove(DIAG)
    unlimited_stack()

    results = {}
    stats = {}
    for n in THREADS:
        got = run(exe, n)
        if got is None:
            return 1
        digests, nkids, exestat = got
        results[n] = (digests, nkids)
        stats[n] = exestat
        print(f"{n:>7} {nkids:>8}  " +
              "  ".join(f"{digests[f][0][:12]}({digests[f][1]:>5}r)" for f in FILES))

    ref_d, ref_k = results[THREADS[0]]
    bad = []
    for n in THREADS[1:]:
        got_d, got_k = results[n]
        if got_k != ref_k:
            bad.append(f"child count {got_k} at {n} threads vs {ref_k} at {THREADS[0]}")
        for f in FILES:
            if got_d[f] != ref_d[f]:
                bad.append(f"{f} differs at {n} threads "
                           f"({got_d[f][1]} records, {got_d[f][0][:12]}) vs {THREADS[0]} "
                           f"({ref_d[f][1]} records, {ref_d[f][0][:12]})")

    if bad:
        for b in bad:
            print(f"  [FAIL] {b}")
        #> Attach evidence before blaming the shed drain: BUGS.md O13 is an unexplained
        #  intermittent truncation that produces exactly this signature on UNMODIFIED code.
        for n in THREADS:
            diag(f"exe at launch, {n} threads: {stats[n]}")
        diag("identical exe across all three launches: "
             f"{len(set(stats.values())) == 1}")
        print("\n[FAIL] KHRT output depends on the thread count. A record-count change points "
              "at the shed drain; equal counts with a different hash points at a TORN record "
              "from the unsynchronised unitTraj/unitExit writes.")
        print(f"  [diag] see {DIAG}; if record counts are a FRACTION of ~10109 this is likely "
              "BUGS.md O13, not a thread bug -- re-run this gate in isolation before believing it.")
        return 1

    print(f"\n[PASS] identical record multisets and child counts across "
          f"OMP_NUM_THREADS = {', '.join(map(str, THREADS))}.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
