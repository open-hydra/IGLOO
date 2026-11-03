#!/usr/bin/env python3
"""psd-bare-token: bc.txt column 8 names a size-distribution file whose name holds neither '.' nor '/'.

A column-8 token that is not a law name is a size-distribution file when one exists under it (as
written, relative to the case directory, then INPUT/ + its basename), whatever its shape. Here the
token is `psd`, a file in the case directory holding one row (radius 12 um): a one-row table is a
Dirac at twice its radius and draws nothing, so every parcel is injected at exactly 24 um -- not the
11.89 um of column 6, which a table ignores.

  L  the fixture: column 8 = `psd` on every inlet line, no '.' or '/'; ./psd holds the typed row and
     INPUT/psd does not exist (the case-directory name is the one resolved)
  N  25 parcels (IDs 1..25), each with an exit record
  D  every parcel's first trajectory row (written at injection) carries d = 2.4e-5 m, to the E13.6
     print's half-ULP
  W  one notice that column 6 lies outside the table, naming `psd`
"""
import os
import sys

TRAJ, EXIT, BC, LOG, TABLE = "OUTPUT/trajectories-A.dat", "OUTPUT/outloc-A.dat", "INPUT/bc.txt", "run_out.txt", "psd"
R1 = 1.2e-5               # the table's single radius [m]
D_INJ = 2.0 * R1          # the injected diameter [m]
RP_COL6 = 5.945e-6        # bc.txt column 6, ignored by a table
NPARCEL = 25
EPS_R = 0.5e-6            # E13.6 relative half-ULP


def rows(path):
    for line in open(path):
        try:
            yield [float(v) for v in line.split()]
        except ValueError:
            continue


def main():
    for f in (TRAJ, EXIT, BC, LOG, TABLE):
        if not os.path.exists(f):
            print(f"[FAIL] {f} missing")
            return 1
    log = open(LOG, errors="replace").read()
    lines = open(BC).read().splitlines()
    payload = [l.split() for l in lines if len(l.split()) == 9]
    tokens = {p[7] for p in payload}
    table = [tuple(float(v) for v in l.split()) for l in open(TABLE)
             if l.strip() and l.strip()[0] not in "#!;"]

    first = {}
    for r in rows(TRAJ):
        if len(r) >= 10:
            pid = int(r[9])
            if pid not in first or r[0] < first[pid][0]:
                first[pid] = r
    exits = {int(r[8]) for r in rows(EXIT) if len(r) >= 9}
    ids = sorted(first)
    worst = max((abs(first[i][7] - D_INJ) / D_INJ for i in ids), default=float("inf"))

    checks = [
        (f"L  column 8 = 'psd' on all {len(payload)} inlet lines (no '.' or '/'), column 6 = {RP_COL6:.5E}",
         tokens == {"psd"} and len(payload) == NPARCEL and all(float(p[5]) == RP_COL6 for p in payload)),
        ("L  ./psd holds the one typed row (radius 1.2e-5, fraction 1.0); INPUT/psd does not exist",
         table == [(R1, 1.0)] and not os.path.exists("INPUT/psd")),
        (f"N  {len(ids)} parcels, IDs 1..{NPARCEL}, each with an exit record",
         ids == list(range(1, NPARCEL + 1)) and exits == set(ids)),
        (f"D  injection diameter 2*r1 = {D_INJ:.6e} m on every first row: worst relative {worst:.2e} <= {EPS_R:.0e}",
         worst <= EPS_R),
        ("W  one notice that column 6 lies outside the table `psd`",
         log.count("[NOTICE]") == 1 and "lies outside psd" in log),
    ]
    fail = 0
    for name, ok in checks:
        print(f"[{'PASS' if ok else 'FAIL'}] {name}")
        fail += (not ok)
    return 1 if fail else 0


if __name__ == "__main__":
    sys.exit(main())
