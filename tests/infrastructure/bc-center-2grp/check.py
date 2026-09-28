#!/usr/bin/env python3
"""Regression gate for the bc_center pinning path with two particle groups (BUGS.md A24).

This case exists to execute `pin_particles_bc_center`, which no other case reaches:
20 e2e cases set `ds` (=> bc_ds), 3 inject from an assigned/DB stream, and every other
`phase.txt` declares one group. The routine held two implicit-SAVE counters, and both
halves of the resulting defect are triggered here (see input.ini for why `fsample = 2`
and `A 2` are both required).

Expectations are built from KNOWN INPUTS, never read back from production output:
the box's face-1 inlet is 5x5 = 25 cells and `fsample = 2` takes every 2nd one, so each
of the 2 groups holds floor(25/2) = 12 particles, 24 in total.

Gates, with the signature each one caught pre-fix (measured 2026-08-07 at d857aba^,
release build, rc=0 and EMPTY STDERR in every case -- none of this announces itself):

  1. population  -- both groups hold exactly 12 particles.
     Pre-fix: group 1 = 12, group 2 = 13. The counting pass carried `nn` across the
     call, so mod(nn,fsample) shifted phase and sized group 2's array from a different
     cell subset than the fill pass walks.
  2. identity    -- each group's ID set is exactly {1..12}.
     Pre-fix: group 1 = {1..12}; group 2 = eleven 0s and a stray 13. `p` resumed at 12,
     so the writes landed at indices 13..24 of a 1:13 array -- one in bounds, ELEVEN
     PAST THE END -- and group 2's real IDs were never assigned, keeping their default 0.
     NOTE: gate on this, not on "IDs are 1..N". A record-count-only gate is vacuous at
     fsample = 1, where both groups return 25 either way.
  3. clone       -- group 2's exit records equal group 1's, particle for particle.
     The two groups see identical cells, properties and physics, so group 2 must be a
     faithful clone. This is the strongest available statement that group 2 was pinned
     correctly rather than merely counted correctly.
  4. integration -- every particle actually flies: >= MIN_TRAJ trajectory records and an
     exit on the outflow plane. A pinning bug that injects everything out of domain
     still exits 0, so "the solver ran" is never the criterion.
"""
import sys
from collections import defaultdict

TRAJ   = "OUTPUT/trajectories-A.dat"
OUTLOC = "OUTPUT/outloc-A.dat"

# ---- known test inputs (NOT read from production output) ---------------------
N_GROUPS   = 2      # INPUT/phase.txt: "A 2"
N_INLET    = 25     # face-1 inlet of the box mesh: 5 x 5 cells
FSAMPLE    = 2      # [IGLOO-BC] fsample
N_PER_GRP  = N_INLET // FSAMPLE     # = 12

X_OUT      = 0.15   # outflow plane (box x-extent)
TOL_X      = 1.0e-6
TOL_CLONE  = 1.0e-9 # groups are bit-identical in practice; slack absorbs a formatting ULP
MIN_TRAJ   = 5      # a particle that dies at injection writes ~1 record


def load_zoned(path, ncols):
    """Return {zone_index: [row, ...]} for a Tecplot point file with one zone per group."""
    zones = defaultdict(list)
    z = 0
    with open(path) as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith("variables"):
                continue
            if s.startswith("Zone"):
                z += 1
                zones[z] = []
                continue
            c = s.split()
            if len(c) != ncols or z == 0:
                continue
            try:
                zones[z].append([float(v) for v in c[:-1]] + [int(c[-1])])
            except ValueError:
                continue
    return zones


def check_population(exits):
    """Gate 1: group count and per-group population."""
    ok = len(exits) == N_GROUPS
    print(f"groups pinned: {len(exits)} (need {N_GROUPS})  [{'PASS' if ok else 'FAIL'}]")
    for g in sorted(exits):
        n = len(exits[g])
        good = n == N_PER_GRP
        ok &= good
        print(f"  group {g}: {n} particles (need {N_PER_GRP} = {N_INLET}//{FSAMPLE})"
              f"  [{'PASS' if good else 'FAIL'}]")
    if len(exits) == N_GROUPS and len({len(v) for v in exits.values()}) != 1:
        print("  [FAIL] groups disagree on population -- the counting pass carried state")
        ok = False
    return 0 if ok else 1


def check_identity(exits):
    """Gate 2: each group's IDs are exactly 1..N_PER_GRP."""
    want = set(range(1, N_PER_GRP + 1))
    ok = True
    for g in sorted(exits):
        ids = [r[-1] for r in exits[g]]
        got = set(ids)
        good = got == want and len(ids) == len(got)
        ok &= good
        if good:
            print(f"  group {g}: IDs 1..{N_PER_GRP}, no duplicates  [PASS]")
        else:
            nz = sum(1 for i in ids if i == 0)
            extra = f"  ({nz} of {len(ids)} are ID 0 -- writes landed out of bounds)" if nz else ""
            print(f"  group {g}: IDs {sorted(got)}{extra}  [FAIL]")
    return 0 if ok else 1


def check_clone(exits):
    """Gate 3: group 2 is a faithful clone of group 1 (same cells => same exits)."""
    if len(exits) < 2:
        print("clone: skipped -- fewer than 2 groups  [FAIL]")
        return 1
    ref = sorted(exits[1], key=lambda r: (r[1], r[2]))   # key on (y, z) injection cell
    ok = True
    for g in sorted(exits):
        if g == 1:
            continue
        rows = sorted(exits[g], key=lambda r: (r[1], r[2]))
        if len(rows) != len(ref):
            print(f"clone: group {g} has {len(rows)} rows vs group 1's {len(ref)}  [FAIL]")
            ok = False
            continue
        worst, where = 0.0, None
        for a, b in zip(ref, rows):
            for k in range(len(a) - 1):        # every column except ID
                d = abs(a[k] - b[k])
                if d > worst:
                    worst, where = d, f"row y={a[1]:.6f} z={a[2]:.6f} col {k}"
        good = worst <= TOL_CLONE
        ok &= good
        print(f"clone: group {g} vs group 1, max |delta| = {worst:.3e} (tol {TOL_CLONE:.0e})"
              f"  [{'PASS' if good else 'FAIL'}]  {where if not good else ''}")
    return 0 if ok else 1


def check_integration(exits, traj):
    """Gate 4: the particles fly and leave through the outflow plane."""
    ok = True
    for g in sorted(exits):
        counts = defaultdict(int)
        for r in traj.get(g, []):
            counts[r[-1]] += 1
        dead = [i for i in (r[-1] for r in exits[g]) if counts[i] < MIN_TRAJ]
        stuck = [r[-1] for r in exits[g] if abs(r[0] - X_OUT) > TOL_X]
        nmin = min(counts.values()) if counts else 0
        good = not dead and not stuck
        ok &= good
        print(f"  group {g}: min {nmin} trajectory records/particle (need >= {MIN_TRAJ}), "
              f"{len(exits[g]) - len(stuck)}/{len(exits[g])} exit at x={X_OUT}"
              f"  [{'PASS' if good else 'FAIL'}]")
        if dead:
            print(f"    [FAIL] IDs with too few records (died at injection?): {sorted(dead)}")
        if stuck:
            print(f"    [FAIL] IDs not on the outflow plane: {sorted(stuck)}")
    return 0 if ok else 1


def main():
    try:
        exits = load_zoned(OUTLOC, 9)
        traj = load_zoned(TRAJ, 10)
    except FileNotFoundError as e:
        print(f"[FAIL] {e.filename} not found -- did the solver run?")
        return 1
    if not exits:
        print(f"[FAIL] no exit records in {OUTLOC}")
        return 1

    print(f"oracle: {N_GROUPS} groups x {N_INLET}//{FSAMPLE} = {N_PER_GRP} particles "
          f"= {N_GROUPS * N_PER_GRP} total, from known inputs\n")
    print("population:")
    rc = check_population(exits)
    print("identity:")
    rc |= check_identity(exits)
    rc |= check_clone(exits)
    print("integration:")
    rc |= check_integration(exits, traj)

    print("\n[PASS] bc_center pinned both groups independently and correctly." if rc == 0
          else "\n[FAIL] bc_center multi-group pinning violation(s) -- see above.")
    return rc


if __name__ == "__main__":
    sys.exit(main())
