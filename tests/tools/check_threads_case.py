#!/usr/bin/env python3
"""Thread-count invariance of a serial e2e case: the same run at several OpenMP thread counts.

    check_threads_case.py <case dir> <IGLOO executable> [threads ...]     (default 1 2 5)

Each thread count runs the case's input.ini in its own directory <case>/threads-omp<N>/ (INPUT/ and input.ini
symlinked to the case's, stdout in run_out.txt, stderr in run_err.txt), so it never touches the case's own
OUTPUT/. Gates:
  1. WITNESS -- at one thread the run prints no `OpenMP threads =` line (the banner appears only when n > 1), at
     n > 1 it prints that n: a binary built without OpenMP, or an environment that never reached the solver,
     cannot pass as a thread comparison.
  2. FILES -- the same OUTPUT/*.dat and *.tec at every thread count, every number finite (check_finite.py).
  3. PARTICLE STREAMS -- every .dat as an exact sorted multiset of its lines: the record order belongs to the
     schedule, the records do not (a torn or a schedule-dependent record is a different line).
  4. GRID FIELDS -- every .tec through compare_tec.py at 1e-12 of the field scale (the `!$OMP ATOMIC` sums
     re-associate with the schedule; their floor is <= 4e-15).
  5. ORACLE -- the case's own check.py on the one-thread run.
PASS stdout prints no .tec residual, so it is the same bytes on every passing run.
"""
import glob
import hashlib
import os
import re
import resource
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from compare_tec import compare_tec, DEFAULT_TOL  # noqa: E402

WITNESS = re.compile(r"OpenMP threads =\s*(\d+)")


def unlimited_stack():
    """The e2e helper runs the solver under `ulimit -s unlimited`; so does this gate (best effort)."""
    try:
        soft, hard = resource.getrlimit(resource.RLIMIT_STACK)
        if soft != resource.RLIM_INFINITY:
            resource.setrlimit(resource.RLIMIT_STACK, (hard, hard))
    except (ValueError, OSError):
        pass


def multiset(path):
    """(md5 of the sorted lines, line count)."""
    with open(path, "rb") as fh:
        lines = sorted(fh.read().splitlines())
    h = hashlib.md5()
    for ln in lines:
        h.update(ln + b"\n")
    return h.hexdigest(), len(lines)


def run(case, exe, n):
    """Run the case at n threads in <case>/threads-omp<n>/; the run directory, or None after printing why."""
    d = os.path.join(case, "threads-omp%d" % n)
    shutil.rmtree(d, ignore_errors=True)
    os.makedirs(os.path.join(d, "OUTPUT"))
    os.symlink(os.path.join("..", "INPUT"), os.path.join(d, "INPUT"))
    os.symlink(os.path.join("..", "input.ini"), os.path.join(d, "input.ini"))
    env = dict(os.environ, OMP_NUM_THREADS=str(n), KMP_STACKSIZE="100M", IGLOO_EXE=exe)
    with open(os.path.join(d, "run_out.txt"), "w") as out, open(os.path.join(d, "run_err.txt"), "w") as err:
        rc = subprocess.call([exe], stdout=out, stderr=err, env=env, cwd=d)
    if rc != 0:
        print("[FAIL] solver exit %d at %d thread(s)" % (rc, n))
        return None
    seen = WITNESS.findall(open(os.path.join(d, "run_out.txt"), errors="replace").read())
    if n == 1 and seen:
        print("[FAIL] 1 thread: the run printed 'OpenMP threads = %s'" % seen[0])
        return None
    if n > 1 and seen != [str(n)]:
        print("[FAIL] %d threads: the run printed %s, not 'OpenMP threads = %d' -- not an OpenMP comparison"
              % (n, seen or "no thread line", n))
        return None
    fin = subprocess.call([sys.executable, "-B", os.path.join(HERE, "check_finite.py")], cwd=d)
    if fin != 0:
        print("[FAIL] %d thread(s): a non-finite number in OUTPUT/" % n)
        return None
    return d


def main():
    if len(sys.argv) < 3:
        print("[FAIL] usage: check_threads_case.py <case dir> <IGLOO executable> [threads ...]")
        return 1
    case, exe = os.path.abspath(sys.argv[1]), os.path.abspath(sys.argv[2])
    threads = [int(t) for t in sys.argv[3:]] or [1, 2, 5]
    if threads[0] != 1:
        print("[FAIL] the first thread count must be 1 (the case's check.py runs on it)")
        return 1
    unlimited_stack()
    name = os.path.basename(case)
    print("thread-count gate on %s: OMP_NUM_THREADS = %s" % (name, ", ".join(map(str, threads))))

    dirs, dats, tecs = {}, {}, {}
    for n in threads:
        d = run(case, exe, n)
        if d is None:
            return 1
        dirs[n] = d
        dats[n] = {os.path.basename(p): multiset(p) for p in sorted(glob.glob(os.path.join(d, "OUTPUT", "*.dat")))}
        tecs[n] = sorted(os.path.basename(p) for p in glob.glob(os.path.join(d, "OUTPUT", "*.tec")))
        print("  %d thread(s): %s; %d .tec" % (n, "  ".join("%s=%s(%dL)" % (f, h[:12], nl)
                                                            for f, (h, nl) in sorted(dats[n].items())), len(tecs[n])))
    ref = threads[0]
    bad = []
    if not dats[ref]:
        bad.append("no .dat output at 1 thread: nothing to compare")
    for n in threads[1:]:
        if sorted(dats[n]) != sorted(dats[ref]) or tecs[n] != tecs[ref]:
            bad.append("the files at %d threads differ from those at 1: %s / %s against %s / %s"
                       % (n, sorted(dats[n]), tecs[n], sorted(dats[ref]), tecs[ref]))
            continue
        for f in sorted(dats[ref]):
            if dats[n][f] != dats[ref][f]:
                bad.append("%s: record multiset at %d threads (%d lines, %s) differs from 1 thread (%d lines, %s)"
                           % (f, n, dats[n][f][1], dats[n][f][0][:12], dats[ref][f][1], dats[ref][f][0][:12]))
        for f in tecs[ref]:
            ok, worst, ndiff, nval = compare_tec(os.path.join(dirs[ref], "OUTPUT", f), os.path.join(dirs[n], "OUTPUT", f))
            if not ok:
                bad.append("%s: scale-relative difference %.3e > %.0e at %d threads (%d of %d values)"
                           % (f, worst, DEFAULT_TOL, n, ndiff, nval))
    if bad:
        for b in bad:
            print("  [FAIL] %s" % b)
        print("[FAIL] %s depends on the thread count" % name)
        return 1

    res = subprocess.run([sys.executable, "-B", os.path.join(case, "check.py")], cwd=dirs[ref],
                         env=dict(os.environ, IGLOO_EXE=exe), stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                         universal_newlines=True)
    with open(os.path.join(dirs[ref], "check_out.txt"), "w") as fh:
        fh.write(res.stdout)
    last = (res.stdout.strip().splitlines() or ["(no output)"])[-1]
    if res.returncode != 0:
        print(res.stdout[-2000:])
        print("[FAIL] %s's check.py on the 1-thread run: exit %d" % (name, res.returncode))
        return 1
    print("  check.py on the 1-thread run: %s" % last.strip())
    print("[PASS] %s: %d .dat multisets exact and %d .tec within %.0e of the field scale at %s threads"
          % (name, len(dats[ref]), len(tecs[ref]), DEFAULT_TOL, ", ".join(map(str, threads))))
    return 0


if __name__ == "__main__":
    sys.exit(main())
