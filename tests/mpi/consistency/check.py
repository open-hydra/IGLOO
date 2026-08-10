#!/usr/bin/env python3
"""MPI rank-count consistency gate: one rank and four ranks must produce the same run.

    check.py <mpiexec> <numproc-flag> <path-to-IGLOO>

The other cases under `tests/mpi/` prove that each case's own oracle still passes under mpiexec.
That is necessary and not sufficient: an oracle checks PHYSICS, and the things MPI can break here are
mostly not physics. This gate compares a 1-rank run against a 4-rank run of the same binary, with the
rank count as the only variable, on four axes:

  1. WITNESS -- the 4-rank run must report `MPI ranks = 4` and the 1-rank run must NOT print the line
     at all. Without this the gate is vacuous in the most likely wrong configuration: `bin/IGLOO` is
     ONE link target shared by every build tree, so a USE_MPI=OFF binary can end up here, and under
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

Plus an integration criterion: two runs that both produced nothing would agree perfectly. Both runs
must show particles that actually moved.

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

RANKS = (1, 4)
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


def moved(recs):
    """Max spread in x over all records -- proves the population is not frozen at injection."""
    if not recs:
        return 0.0
    xs = [r[0] for r in recs]
    return max(xs) - min(xs)


def run(mpiexec, npflag, exe, n, keep):
    shutil.rmtree("OUTPUT", ignore_errors=True)
    os.makedirs("OUTPUT", exist_ok=True)
    log = f"run_n{n}.txt"
    with open(log, "w") as out, open(f"err_n{n}.txt", "w") as err:
        rc = subprocess.call([mpiexec, npflag, str(n), exe], stdout=out, stderr=err)
    if rc != 0:
        print(f"[FAIL] mpiexec exit {rc} at n={n}")
        print(open(log, errors='replace').read()[-800:])
        return None
    text = open(log, errors='replace').read()

    #> Axis 1: the witness. Absent at n=1 by design (the banner prints only when decomposed).
    banner = f"MPI ranks = {n}"
    if n > 1 and banner not in text:
        print(f"[FAIL] n={n}: solver never reported '{banner}'. Either bin/IGLOO is not an MPI "
              f"build (it is one link target shared by every build tree) or MPI_COMM_SIZE "
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
        snap[os.path.basename(p)] = (multiset_hash(recs), len(recs), nline, nzone, moved(recs))
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

    print(f"MPI consistency gate: n = {' vs '.join(str(n) for n in RANKS)}, "
          f"rank count is the only variable")
    snaps, keeps = {}, {}
    for n in RANKS:
        keeps[n] = f"tec_n{n}"
        shutil.rmtree(keeps[n], ignore_errors=True)
        os.makedirs(keeps[n], exist_ok=True)
        got = run(mpiexec, npflag, exe, n, keeps[n])
        if got is None:
            return 1
        snaps[n] = got
        print(f"  n={n}: " + "  ".join(
            f"{name}={h}({nr}r,{nl}L,{nz}Z)" for name, (h, nr, nl, nz, _) in sorted(got.items())))

    ref = RANKS[0]
    bad = []

    #> Integration criterion: two empty runs would agree perfectly.
    for n in RANKS:
        tot = sum(nr for _, nr, _, _, _ in snaps[n].values())
        if tot < MIN_RECORDS:
            bad.append(f"n={n} produced only {tot} records (< {MIN_RECORDS}) -- too thin to compare")
        if not any(mv > MIN_MOVE for *_, mv in snaps[n].values()):
            bad.append(f"n={n}: no particle moved more than {MIN_MOVE} m -- the population never "
                       f"integrated, so any agreement here is vacuous")

    #> Axes 2 (second half) and 3.
    missing = set(snaps[ref]) ^ set(snaps[RANKS[-1]])
    for name in sorted(missing):
        bad.append(f"{name}: present at one rank count and not the other")
    for n in RANKS[1:]:
        for name, (h, nr, nl, nz, _) in sorted(snaps[n].items()):
            if name not in snaps[ref]:
                continue
            rh, rnr, rnl, rnz, _ = snaps[ref][name]
            if nl != rnl:
                bad.append(f"{name}: {nl} lines at n={n} vs {rnl} at n={ref} "
                           f"(a duplicated or dropped header/zone)")
            if nz != rnz:
                bad.append(f"{name}: {nz} Zone records at n={n} vs {rnz} at n={ref}")
            if h != rh:
                bad.append(f"{name}: record multiset differs -- {nr} records {h} at n={n} "
                           f"vs {rnr} records {rh} at n={ref}")

    #> Axis 4.
    ntec = 0
    for src in sorted(glob.glob(os.path.join(keeps[ref], "*.tec"))):
        name = os.path.basename(src)
        for n in RANKS[1:]:
            other = os.path.join(keeps[n], name)
            if not os.path.exists(other):
                bad.append(f"{name}: written at n={ref} but not at n={n}")
                continue
            ok, worst, ndiff, nval = compare_tec(src, other)
            ntec += 1
            if not ok:
                bad.append(f"{name}: scale-relative error {worst:.3e} exceeds {DEFAULT_TOL:.0e} "
                           f"between n={ref} and n={n} ({ndiff} of {nval} values differ)")
    if ntec == 0:
        bad.append("no .tec field was compared -- this fixture is supposed to write source/euler "
                   "output, so axis 4 silently covered nothing")

    if bad:
        for b in bad:
            print(f"  [FAIL] {b}")
        return 1

    print(f"  [PASS] {len(snaps[ref])} .dat multisets exact, layout identical, "
          f"{ntec} .tec field(s) within scale-relative {DEFAULT_TOL:.0e}, across n = "
          f"{', '.join(str(n) for n in RANKS)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
