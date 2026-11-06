#!/usr/bin/env python3
"""Cross-rank reduction determinism: at a fixed rank count the grid fields are the same BYTES whatever
collective algorithm the MPI library picks, and stay within the tolerance of the one-rank run.

    check.py <mpiexec> <numproc-flag> <path-to-IGLOO>

The partial source/euler sums of the ranks are added in rank order, ((v0 + v1) + v2) + ..., and the
total is broadcast, so the association never depends on the library. A library reduction does: the
library chooses its own association per message size and rank count (Intel MPI exposes the choice as
I_MPI_ADJUST_ALLREDUCE, Open MPI as coll_tuned_allreduce_algorithm), and floating-point addition is not
associative. The fixture is the 1672-parcel
vortex cloud (standard/vortex-cloud-st1): every cell of the cloud's path is fed by parcels of every
rank, which is what lets a different association show in the printed digits (a box fixture with one
parcel per cell row cannot see it). Every run is at one thread per rank: the !$OMP ATOMIC deposits
inside a rank are themselves order-dependent with more than one thread.

Per rank count n in (3, 4): run A with the library's default algorithms, run B with its ring allreduce:
I_MPI_ADJUST_ALLREDUCE=8 for Intel MPI, OMPI_MCA_coll_tuned_use_dynamic_rules=1 with
OMPI_MCA_coll_tuned_allreduce_algorithm=4 for Open MPI (each library ignores the other's variables). The ring
associates differently from Intel MPI's default on the reference host, and from Open MPI 4.1's at n = 4
(Rabenseifner for this fixture's message sizes; at n = 3 Open MPI's default is the ring itself).
Gates: rank witness; A and B .tec byte-identical; the .dat files of A, B and n = 1 byte-identical; A within
scale-relative 1e-12 of n = 1. PASS stdout carries no residuals.
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

RANKS = (3, 4)
VARIANT = {'I_MPI_ADJUST_ALLREDUCE': '8',                     # Intel MPI: ring
           'OMPI_MCA_coll_tuned_use_dynamic_rules': '1',      # Open MPI: obey the forced algorithm below
           'OMPI_MCA_coll_tuned_allreduce_algorithm': '4'}    # Open MPI: ring


def md5(path):
    return hashlib.md5(open(path, 'rb').read()).hexdigest()


def records(path):
    out = []
    with open(path, errors='replace') as fh:
        for line in fh:
            s = line.strip()
            if not s or s.lower().startswith(('variables', 'zone')):
                continue
            try:
                out.append(tuple(float(t) for t in s.split()))
            except ValueError:
                continue
    return collections.Counter(out)


def run(mpiexec, npflag, exe, n, tag, extra_env):
    shutil.rmtree('OUTPUT', ignore_errors=True)
    os.makedirs('OUTPUT')
    env = dict(os.environ, OMP_NUM_THREADS='1', KMP_STACKSIZE='100M')
    for k in VARIANT:
        env.pop(k, None)
    env.update(extra_env)
    log = f'run_{tag}.txt'
    with open(log, 'w') as out, open(f'err_{tag}.txt', 'w') as err:
        rc = subprocess.call([mpiexec, npflag, str(n), exe], stdout=out, stderr=err, env=env)
    text = open(log, errors='replace').read()
    if rc != 0:
        print(f'[FAIL] {tag}: mpiexec exit {rc}')
        print(text[-800:])
        return None
    if n > 1 and f'MPI ranks = {n}' not in text:
        print(f"[FAIL] {tag}: solver never reported 'MPI ranks = {n}' -- not an MPI build?")
        return None
    keep = f'tec_{tag}'
    shutil.rmtree(keep, ignore_errors=True)
    shutil.copytree('OUTPUT', keep)
    tec = sorted(os.path.basename(p) for p in glob.glob(os.path.join(keep, '*.tec')))
    dat = sorted(os.path.basename(p) for p in glob.glob(os.path.join(keep, '*.dat')))
    if not tec or not dat:
        print(f'[FAIL] {tag}: no .tec or no .dat output')
        return None
    return keep, tec, dat


def main():
    if len(sys.argv) < 4:
        print('[FAIL] usage: check.py <mpiexec> <numproc-flag> <path-to-IGLOO>')
        return 1
    mpiexec, npflag, exe = sys.argv[1:4]
    print(f'MPI reduction determinism: n = {", ".join(str(n) for n in RANKS)}, default vs '
          f'{" ".join(f"{k}={v}" for k, v in VARIANT.items())}, one thread per rank')
    ref = run(mpiexec, npflag, exe, 1, 'n1', {})
    if ref is None:
        return 1
    bad = []
    for n in RANKS:
        a = run(mpiexec, npflag, exe, n, f'n{n}a', {})
        b = run(mpiexec, npflag, exe, n, f'n{n}b', VARIANT)
        if a is None or b is None:
            return 1
        for name in sorted(set(ref[1]) | set(a[1]) | set(b[1])):
            pa, pb, pr = (os.path.join(x[0], name) for x in (a, b, ref))
            if not (os.path.exists(pa) and os.path.exists(pb) and os.path.exists(pr)):
                bad.append(f'n={n} {name}: missing in one of the runs')
                continue
            if md5(pa) != md5(pb):
                ok, worst, nd, nv = compare_tec(pa, pb)
                bad.append(f'n={n} {name}: bytes differ between the default and the variant algorithm '
                           f'({nd} of {nv} values, scale-relative {worst:.3e})')
            ok, worst, nd, nv = compare_tec(pr, pa)
            if not ok:
                bad.append(f'n={n} {name}: scale-relative {worst:.3e} vs n=1 exceeds {DEFAULT_TOL:.0e}')
        for name in sorted(set(ref[2]) | set(a[2]) | set(b[2])):
            paths = [os.path.join(x[0], name) for x in (ref, a, b)]
            if not all(os.path.exists(p) for p in paths):
                bad.append(f'n={n} {name}: missing in one of the runs')
                continue
            hr, ha, hb = (md5(p) for p in paths)
            if not (hr == ha == hb):
                cr, ca, cb = (records(p) for p in paths)
                bad.append(f'n={n} {name}: bytes differ from n=1 ({sum(cr.values())}, {sum(ca.values())}, '
                           f'{sum(cb.values())} records; multisets {"equal" if cr == ca == cb else "differ"})')
    if bad:
        for m in bad:
            print(f'  [FAIL] {m}')
        return 1
    print(f'  [PASS] {len(ref[1])} .tec field(s) byte-identical across the algorithm choice at '
          f'n = {", ".join(str(n) for n in RANKS)} and within {DEFAULT_TOL:.0e} of n = 1; '
          f'{len(ref[2])} .dat byte-identical to n = 1')
    return 0


if __name__ == '__main__':
    sys.exit(main())
