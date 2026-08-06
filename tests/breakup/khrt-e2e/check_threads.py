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
import hashlib
import os
import shutil
import subprocess
import sys

THREADS = (1, 2, 4)
FILES = ("OUTPUT/trajectories-A.dat", "OUTPUT/outloc-A.dat")


def digest(path):
    """Hash of the record multiset: sorted lines, so ordering is ignored but content is not."""
    with open(path) as fh:
        lines = sorted(fh.readlines())
    h = hashlib.md5()
    for ln in lines:
        h.update(ln.encode())
    return h.hexdigest(), len(lines)


def run(exe, nthreads):
    shutil.rmtree("OUTPUT", ignore_errors=True)
    os.makedirs("OUTPUT", exist_ok=True)
    env = dict(os.environ, OMP_NUM_THREADS=str(nthreads), KMP_STACKSIZE="100M")
    with open("run_out.txt", "w") as out, open("run_err.txt", "w") as err:
        rc = subprocess.call([exe], stdout=out, stderr=err, env=env)
    if rc >= 128:
        print(f"[FAIL] solver died on signal (exit {rc}) at OMP_NUM_THREADS={nthreads}")
        return None
    if rc != 0:
        print(f"[FAIL] solver exit {rc} at OMP_NUM_THREADS={nthreads}")
        return None
    nkids = 0
    for line in open("run_out.txt"):
        if "number of children" in line:
            nkids += int(line.split("=")[-1])
    try:
        return {f: digest(f) for f in FILES}, nkids
    except FileNotFoundError as exc:
        print(f"[FAIL] {exc.filename} missing at OMP_NUM_THREADS={nthreads}")
        return None


def main():
    if len(sys.argv) < 2:
        print("[FAIL] usage: check_threads.py <path-to-IGLOO-binary>")
        return 1
    exe = sys.argv[1]

    print("KHRT thread-invariance gate (sorted-multiset compare):")
    print(f"{'threads':>7} {'children':>8}  " +
          "  ".join(f"{os.path.basename(f):>24}" for f in FILES))

    results = {}
    for n in THREADS:
        got = run(exe, n)
        if got is None:
            return 1
        digests, nkids = got
        results[n] = (digests, nkids)
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
        print("\n[FAIL] KHRT output depends on the thread count. A record-count change points "
              "at the shed drain; equal counts with a different hash points at a TORN record "
              "from the unsynchronised unitTraj/unitExit writes.")
        return 1

    print(f"\n[PASS] identical record multisets and child counts across "
          f"OMP_NUM_THREADS = {', '.join(map(str, THREADS))}.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
