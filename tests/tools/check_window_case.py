#!/usr/bin/env python3
"""The particle files' ordered window under a small memory budget: the same bytes as at one thread.

    check_window_case.py <case dir> <IGLOO executable> [threads] [budget]     (default 5 threads, 65536 bytes)

The case's input.ini, with `timing = on` added to [IGLOO-General], runs three times in its own directories under the
case (INPUT/ symlinked, stdout in run_out.txt), so it never touches the case's own OUTPUT/:
window-omp1/ at one thread, window-default/ at <threads> threads, window-small/ at <threads> threads with
IGLOO_IOBUF_BYTES = <budget> (the window's memory budget; the default is 64 MiB). Gates:
  1. WITNESS -- each run prints its `OpenMP threads =` line (none at one thread) and one timing-report line
     `particle files  window budget B ...` with the budget it ran with.
  2. ORACLE -- at one thread every block arrives in its turn: held max 0, spilled 0.
  3. NON-VACUITY -- under the small budget a block arrived before its turn and did not fit: spilled > 0.
  4. BOUND -- under the small budget the bytes held never exceeded it: held max <= budget.
  5. FILES -- every OUTPUT/*.dat of both <threads>-thread runs byte for byte the one-thread run's.
"""
import glob
import hashlib
import os
import re
import resource
import shutil
import subprocess
import sys

WITNESS = re.compile(r"OpenMP threads =\s*(\d+)")
WINDOW = re.compile(r"\[IGLOO-timing\] particle files  window budget (\d+) B  held max (\d+) B  spilled (\d+) B"
                    r"  largest block (\d+) B")
DEFAULT_BUDGET = 64 * 1024 * 1024


def unlimited_stack():
    """The e2e helper runs the solver under `ulimit -s unlimited`; so does this gate (best effort)."""
    try:
        soft, hard = resource.getrlimit(resource.RLIMIT_STACK)
        if soft != resource.RLIM_INFINITY:
            resource.setrlimit(resource.RLIMIT_STACK, (hard, hard))
    except (ValueError, OSError):
        pass


def run(case, exe, name, n, budget):
    """Run the case in <case>/<name>/ at n threads (budget: IGLOO_IOBUF_BYTES, or None); (dir, window figures) or
    None after printing why."""
    d = os.path.join(case, name)
    shutil.rmtree(d, ignore_errors=True)
    os.makedirs(os.path.join(d, "OUTPUT"))
    os.symlink(os.path.join("..", "INPUT"), os.path.join(d, "INPUT"))
    lines = []
    for ln in open(os.path.join(case, "input.ini")).read().split("\n"):
        lines.append(ln)
        if ln.strip() == "[IGLOO-General]":
            lines.append("timing = on")
    open(os.path.join(d, "input.ini"), "w").write("\n".join(lines))
    env = dict(os.environ, OMP_NUM_THREADS=str(n), KMP_STACKSIZE="100M", IGLOO_EXE=exe)
    env.pop("IGLOO_IOBUF_BYTES", None)
    if budget is not None:
        env["IGLOO_IOBUF_BYTES"] = str(budget)
    with open(os.path.join(d, "run_out.txt"), "w") as out, open(os.path.join(d, "run_err.txt"), "w") as err:
        rc = subprocess.call([exe], stdout=out, stderr=err, env=env, cwd=d)
    if rc != 0:
        print("[FAIL] %s: solver exit %d" % (name, rc))
        return None
    text = open(os.path.join(d, "run_out.txt"), errors="replace").read()
    seen = WITNESS.findall(text)
    if (n == 1 and seen) or (n > 1 and seen != [str(n)]):
        print("[FAIL] %s: the run printed %s at %d thread(s) -- not this thread count"
              % (name, seen or "no thread line", n))
        return None
    win = WINDOW.findall(text)
    if len(win) != 1:
        print("[FAIL] %s: %d particle-files timing lines, expected 1" % (name, len(win)))
        return None
    fig = [int(v) for v in win[0]]
    want = DEFAULT_BUDGET if budget is None else budget
    if fig[0] != want:
        print("[FAIL] %s: window budget %d B, expected %d B" % (name, fig[0], want))
        return None
    return d, fig


def filehash(path):
    """md5 of the file's bytes."""
    with open(path, "rb") as fh:
        return hashlib.md5(fh.read()).hexdigest()


def main():
    if len(sys.argv) < 3:
        print("[FAIL] usage: check_window_case.py <case dir> <IGLOO executable> [threads] [budget]")
        return 1
    case, exe = os.path.abspath(sys.argv[1]), os.path.abspath(sys.argv[2])
    nthr = int(sys.argv[3]) if len(sys.argv) > 3 else 5
    budget = int(sys.argv[4]) if len(sys.argv) > 4 else 65536
    unlimited_stack()
    name = os.path.basename(case)
    print("ordered-window gate on %s: 1 thread, %d threads at the default budget and at %d B" % (name, nthr, budget))

    runs = {}
    for key, n, b in (("window-omp1", 1, None), ("window-default", nthr, None), ("window-small", nthr, budget)):
        r = run(case, exe, key, n, b)
        if r is None:
            return 1
        runs[key] = r
        print("  %s: budget %d B, held max %d B, spilled %d B, largest block %d B" % ((key,) + tuple(r[1])))
    bad = []
    one = runs["window-omp1"][1]
    if one[1] != 0 or one[2] != 0:
        bad.append("one thread held %d B and spilled %d B: a block arrived out of its turn" % (one[1], one[2]))
    small = runs["window-small"][1]
    if small[2] <= 0:
        bad.append("the %d B budget spilled nothing: the spill path was not exercised" % budget)
    if small[1] > budget:
        bad.append("the %d B budget held %d B" % (budget, small[1]))
    ref = {os.path.basename(p): filehash(p)
           for p in sorted(glob.glob(os.path.join(runs["window-omp1"][0], "OUTPUT", "*.dat")))}
    if not ref:
        bad.append("no .dat output at one thread: nothing to compare")
    for key in ("window-default", "window-small"):
        got = {os.path.basename(p): filehash(p)
               for p in sorted(glob.glob(os.path.join(runs[key][0], "OUTPUT", "*.dat")))}
        if sorted(got) != sorted(ref):
            bad.append("%s wrote %s, one thread %s" % (key, sorted(got), sorted(ref)))
            continue
        for f in sorted(ref):
            if got[f] != ref[f]:
                bad.append("%s: %s differs from the one-thread file (%s against %s)" % (key, f, got[f][:12], ref[f][:12]))
    if bad:
        for b in bad:
            print("  [FAIL] %s" % b)
        print("[FAIL] %s: %d check(s) of the ordered window failed" % (name, len(bad)))
        return 1
    print("[PASS] %s: %d .dat byte-identical at 1 and %d threads, the %d B budget spilled and held at most %d B"
          % (name, len(ref), nthr, budget, budget))
    return 0


if __name__ == "__main__":
    sys.exit(main())
