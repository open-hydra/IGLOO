#!/usr/bin/env python3
"""MPI rank-count and thread-count consistency gate: every (ranks, threads) arm must produce the same run.

    check.py <mpiexec> <numproc-flag> <path-to-IGLOO>

The other cases under `tests/mpi/` prove that each case's own oracle still passes under mpiexec.
That is necessary and not sufficient: an oracle checks PHYSICS, and the things MPI can break here are
mostly not physics. This gate compares runs of the same binary over a rank axis (1 vs 4 ranks at the
ctest thread count) and a thread axis (4 ranks at 1 and at 5 threads per rank) against the 1-rank
run, on four axes:

  1. WITNESS -- the 4-rank run must report `MPI ranks = 4` and the 1-rank run must NOT print the line
     at all. Without this the gate is vacuous in the most likely wrong configuration: a USE_MPI=OFF
     binary passed as the executable (a serial tree's, or bin/IGLOO from the last tree built), and under
     `mpiexec -n 4` it runs four independent full sweeps that clobber one another's OUTPUT/. Every
     after-the-fact file check then passes while nothing was decomposed.
  2. LAYOUT -- no `*.rank<r>.dat` shard survives (the in-code merge ran and deleted them), and the
     merged files carry the same total line count and the same number of `Zone` records as the
     1-rank run. A content-only check does NOT cover this: duplicating the zone header per rank left
     a case's own check.py passing (measured 2026-08-10, Phase 4).
  3. PARTICLE STREAMS -- `.dat` files compare as an exact sorted MULTISET. Record order is
     nondeterministic by design (OMP writes as particles finish, and the merge emits rank-blocked
     data), but the set of records is exact: no tolerance.
  4. GRID FIELDS -- `.tec` files compare on SCALE-relative error (tools/compare_tec.py). These are
     accumulated under `!$OMP ATOMIC UPDATE`, so changing the rank count re-partitions a
     non-associative sum; that is the same FP-reassociation class already floored at <= 4e-15
     run-to-run.

Plus an integration criterion: two runs that both produced nothing would agree perfectly. Every run
must show particles that actually moved: the records of all its .dat files together span more than
MIN_MOVE in x (on a fixture one cell long a parcel's injection row and its exit lie in different files).

The thread axis exists because rank invariance holds "by construction" (per-parcel RNG streams), and a
sampler that drew from a shared stream inside the parallel region would break the thread count first;
at a fixed thread count the rank axis cannot see that.

PASS stdout is deliberately reproducible -- no `.tec` residuals or diff counts, which vary with
thread interleaving (BUGS.md O2: an irreproducible gate stdout pollutes every future byte-inert A/B).
Fully verbose on FAIL, which is when the numbers are wanted.
"""
import collections
import glob
import hashlib
import os
import shutil
import subprocess
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'tools'))
from compare_tec import compare_tec, DEFAULT_TOL          # noqa: E402

#> (ranks, threads per rank); None = the ctest environment's OMP_NUM_THREADS. The first arm is the reference.
ARMS = ((1, None), (4, None), (4, 1), (4, 5))
MIN_RECORDS = 20          #> below this the comparison is too thin to mean anything
MIN_MOVE = 1.0e-9         #> metres; a particle that moved less never left its injection point


def data_lines(path):
    """(data records as float tuples, total line count, zone-header count)."""
    recs, nline, nzone = [], 0, 0
    with open(path, errors='replace') as fh:
        for line in fh:
            nline += 1
            s = line.strip()
            if not s:
                continue
            if s.lower().startswith(('variables', 'zone')):
                if s.lower().startswith('zone'):
                    nzone += 1
                continue
            try:
                recs.append(tuple(float(t) for t in s.split()))
            except ValueError:
                continue
    return recs, nline, nzone


def multiset_hash(recs):
    h = hashlib.md5()
    for r in sorted(recs):
        h.update(repr(r).encode())
    return h.hexdigest()[:12]


def x_range(recs):
    """(min x, max x) over the records, None without records: the spread over a run's files proves the
    population was not frozen at injection."""
    if not recs:
        return None
    xs = [r[0] for r in recs]
    return min(xs), max(xs)


def tag(arm):
    n, t = arm
    return f"n{n}" if t is None else f"n{n}_t{t}"


def label(arm):
    n, t = arm
    return f"n={n}" if t is None else f"n={n},omp={t}"


def run(mpiexec, npflag, exe, arm, keep):
    n, t = arm
    shutil.rmtree("OUTPUT", ignore_errors=True)
    os.makedirs("OUTPUT", exist_ok=True)
    env = dict(os.environ)
    if t is not None:
        env["OMP_NUM_THREADS"] = str(t)
    log = f"run_{tag(arm)}.txt"
    with open(log, "w") as out, open(f"err_{tag(arm)}.txt", "w") as err:
        rc = subprocess.call([mpiexec, npflag, str(n), exe], stdout=out, stderr=err, env=env)
    if rc != 0:
        print(f"[FAIL] mpiexec exit {rc} at {label(arm)}")
        print(open(log, errors='replace').read()[-800:])
        return None
    text = open(log, errors='replace').read()

    #> Axis 1: the witness. Absent at n=1 by design (the banner prints only when decomposed).
    banner = f"MPI ranks = {n}"
    if n > 1 and banner not in text:
        print(f"[FAIL] n={n}: solver never reported '{banner}'. Either the executable is not an MPI "
              f"build (a serial tree's) or MPI_COMM_SIZE "
              f"disagrees with mpiexec -- in the first case nothing was decomposed and this "
              f"comparison would have been vacuous.")
        return None
    if n == 1 and "MPI ranks" in text:
        print(f"[FAIL] n=1: solver printed an 'MPI ranks' banner, so the serial path is no longer "
              f"byte-identical to a single-rank run.")
        return None

    #> Axis 2, first half: the merge must have consumed every shard.
    left = sorted(glob.glob("OUTPUT/*.rank*.dat"))
    if left:
        print(f"[FAIL] n={n}: {len(left)} rank shard(s) survived the merge: "
              f"{', '.join(os.path.basename(p) for p in left)}")
        return None

    snap = {}
    for p in sorted(glob.glob("OUTPUT/*.dat")):
        recs, nline, nzone = data_lines(p)
        snap[os.path.basename(p)] = (multiset_hash(recs), len(recs), nline, nzone, x_range(recs))
    for p in sorted(glob.glob("OUTPUT/*.tec")):
        dst = os.path.join(keep, os.path.basename(p))
        shutil.copyfile(p, dst)
    if not snap:
        print(f"[FAIL] n={n}: no .dat output at all")
        return None
    return snap


def main():
    if len(sys.argv) < 4:
        print("[FAIL] usage: check.py <mpiexec> <numproc-flag> <path-to-IGLOO>")
        return 1
    mpiexec, npflag, exe = sys.argv[1], sys.argv[2], sys.argv[3]

    print(f"MPI consistency gate: {' vs '.join(label(a) for a in ARMS)}; the rank count, then the "
          f"thread count, is the only variable")
    snaps, keeps = {}, {}
    for arm in ARMS:
        keeps[arm] = f"tec_{tag(arm)}"
        shutil.rmtree(keeps[arm], ignore_errors=True)
        os.makedirs(keeps[arm], exist_ok=True)
        got = run(mpiexec, npflag, exe, arm, keeps[arm])
        if got is None:
            return 1
        snaps[arm] = got
        print(f"  {label(arm)}: " + "  ".join(
            f"{name}={h}({nr}r,{nl}L,{nz}Z)" for name, (h, nr, nl, nz, _) in sorted(got.items())))

    ref = ARMS[0]
    bad = []

    #> Integration criterion: two empty runs would agree perfectly.
    for arm in ARMS:
        tot = sum(nr for _, nr, _, _, _ in snaps[arm].values())
        if tot < MIN_RECORDS:
            bad.append(f"{label(arm)} produced only {tot} records (< {MIN_RECORDS}) -- too thin to compare")
        spans = [xr for *_, xr in snaps[arm].values() if xr is not None]
        if not spans or max(hi for _, hi in spans) - min(lo for lo, _ in spans) <= MIN_MOVE:
            bad.append(f"{label(arm)}: no particle moved more than {MIN_MOVE} m -- the population never "
                       f"integrated, so any agreement here is vacuous")

    #> Axes 2 (second half) and 3.
    for arm in ARMS[1:]:
        for name in sorted(set(snaps[ref]) ^ set(snaps[arm])):
            bad.append(f"{name}: present at {label(ref)} or at {label(arm)}, not both")
        for name, (h, nr, nl, nz, _) in sorted(snaps[arm].items()):
            if name not in snaps[ref]:
                continue
            rh, rnr, rnl, rnz, _ = snaps[ref][name]
            if nl != rnl:
                bad.append(f"{name}: {nl} lines at {label(arm)} vs {rnl} at {label(ref)} "
                           f"(a duplicated or dropped header/zone)")
            if nz != rnz:
                bad.append(f"{name}: {nz} Zone records at {label(arm)} vs {rnz} at {label(ref)}")
            if h != rh:
                bad.append(f"{name}: record multiset differs -- {nr} records {h} at {label(arm)} "
                           f"vs {rnr} records {rh} at {label(ref)}")

    #> Axis 4.
    ntec = 0
    for src in sorted(glob.glob(os.path.join(keeps[ref], "*.tec"))):
        name = os.path.basename(src)
        for arm in ARMS[1:]:
            other = os.path.join(keeps[arm], name)
            if not os.path.exists(other):
                bad.append(f"{name}: written at {label(ref)} but not at {label(arm)}")
                continue
            ok, worst, ndiff, nval = compare_tec(src, other)
            ntec += 1
            if not ok:
                bad.append(f"{name}: scale-relative error {worst:.3e} exceeds {DEFAULT_TOL:.0e} "
                           f"between {label(ref)} and {label(arm)} ({ndiff} of {nval} values differ)")
    if ntec == 0:
        bad.append("no .tec field was compared -- this fixture is supposed to write source/euler "
                   "output, so axis 4 silently covered nothing")

    if bad:
        for b in bad:
            print(f"  [FAIL] {b}")
        return 1

    print(f"  [PASS] {len(snaps[ref])} .dat multisets exact, layout identical, "
          f"{ntec} .tec comparison(s) within scale-relative {DEFAULT_TOL:.0e}, across "
          f"{', '.join(label(a) for a in ARMS)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
