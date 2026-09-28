#!/usr/bin/env python3
"""Two-sweep repeatability oracle. Run from a case directory after tests/support/twosweep.

    check_twosweep.py steady      sweep 1 must REPRODUCE sweep 0
    check_twosweep.py gas-cycle   sweep 1 (U doubled) must DIFFER from sweep 0,
                                  sweep 2 (original field again) must reproduce it

Why this test exists: the rest of the suite performs exactly one sweep, so it is
structurally blind to state that leaks from one solve() into the next. Everything it
checks -- trajectories, exit locations, source field -- is checked against the SAME
code on the SAME input; the only variable is how many times solve() has run.

Two comparison rules, both forced on us by real nondeterminism rather than chosen:

  *.dat  sorted MULTISET equality. Record order is OMP-nondeterministic (the parallel
         loop writes as particles finish), so byte-identity does not hold even between
         two runs of one sweep. The multiset is exact -- no tolerance.

  *.tec  max SCALE-relative difference below TEC_TOL, i.e. max|a-b| normalised by the
         FIELD's magnitude max|a| -- NOT a per-value |a-b|/|a|, which is ill-conditioned
         here (see tools/compare_tec.py for the measurement). The PASS string below says
         "max|rel|" for stdout byte-stability; the computation is and always was
         scale-relative. These fields are filled by
         !$OMP ATOMIC UPDATE, and FP addition is not associative, so the value depends
         on thread interleaving. Measured 2026-08-07 on this host: the run-to-run floor
         at a FIXED sweep is <= 4e-15, while the F7 accumulator-shape bug this test was
         written to catch showed 5.2e-11. TEC_TOL sits ~4 decades from each.

Plus an integration criterion (skill: igloo-verify-integration): a reset that injects
every particle out of the domain also produces two identical empty sweeps and would pass
the comparison vacuously. So sweep 0 must also show particles that actually moved.
"""
import sys, os, glob, re, collections

TEC_TOL = 1.0e-12
MIN_MOVE = 1.0e-9          # metres; below this a particle did not leave its injection point
MOVED_FRAC = 0.9           # allow a few legitimately-stillborn injections

fails, notes = [], []


def records(path):
    """Numeric data rows, as tuples of floats. Skips 'variables=' / 'Zone' headers."""
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
    return out


def tec_values(path):
    vals = []
    with open(path, errors='replace') as fh:
        for line in fh:
            try:
                vals.append(float(line.strip()))
            except ValueError:
                vals.append(None)
    return vals


def cmp_dat(a, b, label):
    ra, rb = records(a), records(b)
    if collections.Counter(ra) == collections.Counter(rb):
        notes.append(f"  {label:<44} {len(ra):>7} records, multiset identical")
        return True
    ca, cb = collections.Counter(ra), collections.Counter(rb)
    fails.append(f"  {label}: {len(ra)} vs {len(rb)} records; "
                 f"{sum((ca - cb).values())} only in first, {sum((cb - ca).values())} only in second")
    return False


def cmp_tec(a, b, label):
    va, vb = tec_values(a), tec_values(b)
    if len(va) != len(vb):
        fails.append(f"  {label}: length {len(va)} vs {len(vb)}")
        return False
    scale = max((abs(x) for x in va if x is not None), default=1.0) or 1.0
    worst, ndiff = 0.0, 0
    for x, y in zip(va, vb):
        if x is None or y is None or x == y:
            continue
        ndiff += 1
        worst = max(worst, abs(x - y) / scale)
    if worst <= TEC_TOL:
        #> Deliberately NOT printing `worst` or `ndiff` on success. Both are set by !$OMP ATOMIC
        #  accumulation order and so vary run to run (measured: khrt 3.9e-16 vs 3.9e-15, vie-plait
        #  1415 vs 1606 diffs, all far inside tolerance). Printing them makes THIS gate's stdout
        #  irreproducible, which pollutes every future byte-inert A/B -- the same trap as BUGS.md
        #  O2. Stable on PASS, fully verbose on FAIL, which is when the numbers are wanted.
        notes.append(f"  {label:<44} within tolerance (max|rel| <= {TEC_TOL:.0e})")
        return True
    fails.append(f"  {label}: max|rel|={worst:.3e} exceeds {TEC_TOL:.0e} ({ndiff} values differ)")
    return False


def differs(a, b, label):
    """For gas-cycle: the perturbed sweep MUST NOT reproduce the baseline."""
    if collections.Counter(records(a)) != collections.Counter(records(b)):
        notes.append(f"  {label:<44} differs, as required")
        return True
    fails.append(f"  {label}: IDENTICAL to the baseline sweep -- the gas refresh "
                 f"never reached the solver")
    return False


def check_moved(mats):
    """Sweep 0 must show particles that left their injection point."""
    for m in mats:
        traj, out = f'OUTPUT/trajectories-{m}.dat', f'OUTPUT/outloc-{m}.dat'
        exits = records(out)
        if not exits:
            fails.append(f"  outloc-{m}.dat: no exit records -- nothing integrated")
            continue
        if not os.path.exists(traj):
            notes.append(f"  outloc-{m}.dat{'':<32} {len(exits)} exits (out-traj off; move check skipped)")
            continue
        first = {}
        for r in records(traj):                       # trajectories: X,Y,Z,...,ID (col 9)
            if len(r) >= 10:
                first.setdefault(int(r[9]), r[0:3])
        moved = tot = 0
        for r in exits:                               # outloc: X,Y,Z,...,ID (col 8)
            if len(r) < 9:
                continue
            inj = first.get(int(r[8]))
            if inj is None:
                continue
            tot += 1
            if sum((r[i] - inj[i]) ** 2 for i in range(3)) ** 0.5 > MIN_MOVE:
                moved += 1
        if tot and moved / tot < MOVED_FRAC:
            fails.append(f"  {m}: only {moved}/{tot} particles left their injection point "
                         f"-- the reset injected them out of domain")
        else:
            notes.append(f"  material {m:<36} {moved}/{tot} particles advanced, {len(exits)} exits")


mode = sys.argv[1] if len(sys.argv) > 1 else 'steady'
mats = sorted(re.match(r'outloc-(.+)\.dat$', os.path.basename(p)).group(1)
              for p in glob.glob('OUTPUT/outloc-*.dat')
              if '-sweep' not in os.path.basename(p))
if not mats:
    print('[FAIL] no OUTPUT/outloc-*.dat from sweep 0 -- the harness did not run')
    sys.exit(1)

check_moved(mats)

if mode == 'steady':
    pairs = [(f'OUTPUT/{k}-{m}.dat', f'OUTPUT/{k}-{m}-sweep1.dat', f'{k}-{m}: sweep1 vs sweep0')
             for m in mats for k in ('outloc', 'trajectories', 'scatter')]
    for a, b, lbl in pairs:
        if os.path.exists(a) and os.path.exists(b):
            cmp_dat(a, b, lbl)
        elif os.path.exists(a) != os.path.exists(b):
            fails.append(f"  {lbl}: only one of the two sweeps wrote this file")
    for a in sorted(glob.glob('OUTPUT/*.tec')):
        if '-sweep' in a:
            continue
        b = a[:-4] + '-sweep1.tec'
        if os.path.exists(b):
            cmp_tec(a, b, f'{os.path.basename(a)}: sweep1 vs sweep0')

elif mode == 'gas-cycle':
    for m in mats:
        differs(f'OUTPUT/outloc-{m}.dat', f'OUTPUT/outloc-{m}-sweep1.dat',
                f'outloc-{m}: sweep1 (U doubled) vs sweep0')
        for k in ('outloc', 'trajectories', 'scatter'):
            a, b = f'OUTPUT/{k}-{m}.dat', f'OUTPUT/{k}-{m}-sweep2.dat'
            if os.path.exists(a) and os.path.exists(b):
                cmp_dat(a, b, f'{k}-{m}: sweep2 (field restored) vs sweep0')
    for a in sorted(glob.glob('OUTPUT/*.tec')):
        if '-sweep' in a:
            continue
        b = a[:-4] + '-sweep2.tec'
        if os.path.exists(b):
            cmp_tec(a, b, f'{os.path.basename(a)}: sweep2 vs sweep0')
else:
    print(f'[FAIL] unknown mode {mode!r}')
    sys.exit(1)

print(f'--- two-sweep repeatability [{mode}] ---')
for n in notes:
    print(n)
if fails:
    print('[FAIL] sweeps disagree:')
    for f in fails:
        print(f)
    sys.exit(1)
print(f'[PASS] {mode}: {len(notes)} comparisons, all within the measured noise floor')
