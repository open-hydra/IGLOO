#!/usr/bin/env python3
"""tc-box-ord2-row: the ord2 source reduction must CONSERVE, end to end.

Telescoping audit. Every kilogram a parcel loses inside the domain must appear in
OUTPUT/source.tec, whatever mesh row the parcel happened to be in:

    sum over cells of wdot  ==  sum over parcels of npdot * (m_inject - m_exit)

Before O28 the ord2 reduction volume-AVERAGED the dual accumulators onto the geo mesh, so a
dual cell was handed on with total weight sum(subVol)/V_geo -- 1 in the interior but 1/2 on a
face row, 1/4 on an edge, 1/8 on a corner. These four parcels ride the j = 1 dual row for
their whole path, so essentially the whole deposit was halved: resid ~ 0.50 against this
oracle. After the fix the weights are a partition of unity and the identity is exact up to
the print truncation.

The oracle reads only trajectory rows and the injected mass flow; it never reads the
Eulerian field it is checking.
"""
import importlib.util
import math
import os
import re
import sys

HERE   = os.path.dirname(os.path.abspath(__file__))
TRAJ   = "OUTPUT/trajectories-A.dat"
OUTLOC = "OUTPUT/outloc-A.dat"
SRC    = "OUTPUT/source.tec"
LOG    = "run_out.txt"

NPART = 4
LX    = 0.15          # box length; the parcels exit through face 2
DX    = 0.0025        # geo cell size in x
#> Measured: GREEN resid 3.779e-06, identical at OMP 1/2/5 (the four parcels ride separate
#> (j,k) columns, so no two threads accumulate into one cell); RED on the pre-O28 binary
#> 5.015e-01. This bound sits 265x above the floor and 500x below the signature. The floor
#> is the E13.6 mass print plus the last-segment closure, both second order on the wet-bulb
#> plateau -- do not loosen past ~1e-2, where the gate stops separating half from whole.
TELESCOPE_TOL = 1.0e-3


def load_trajectories(path):
    """tc-box's own reader, loaded BY PATH: `import check` from a file named check.py
    imports this file itself."""
    spec = importlib.util.spec_from_file_location(
        "tcbox_check", os.path.join(HERE, "..", "tc-box", "check.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod.load_trajectories(path)


def read_source_total(path):
    """Sum wdot over every cell. source.tec is geo-shaped even under ord2 (finalizeSRC
    move_allocs the accumulators down before the write), BLOCK packed, wdot first."""
    lines = open(path).read().splitlines()
    zone = next(l for l in lines if l.strip().startswith("ZONE"))
    I, J, K = (int(re.search(rf"{a}=(\d+)", zone).group(1)) for a in "IJK")
    nn, nc = I*J*K, (I-1)*(J-1)*(K-1)
    vals = " ".join(lines[lines.index(zone)+1:]).split()
    return sum(float(v) for v in vals[3*nn:3*nn + nc]), (I-1, J-1, K-1)


def main():
    for f in (TRAJ, OUTLOC, SRC, LOG):
        if not os.path.exists(f):
            print(f"[FAIL] {f} missing")
            return 1
    log = open(LOG, errors="replace").read()

    parts = load_trajectories(TRAJ)
    mdot_inj = {}
    for line in open(OUTLOC).read().splitlines()[2:]:
        c = line.split()
        if len(c) == 9:
            mdot_inj[int(c[8])] = float(c[6])

    src_total, (nx, ny, nz) = read_source_total(SRC)

    evap = 0.0
    tail_ok = True
    worst_gap = 0.0
    for pid, rows in parts.items():
        rows = sorted(rows, key=lambda r: r[0])
        m_first, m_last = rows[0][8], rows[-1][8]
        x_last = rows[-1][0]
        if pid not in mdot_inj or m_first <= 0.0:
            continue
        npdot = mdot_inj[pid] / m_first
        gap = LX - x_last
        worst_gap = max(worst_gap, gap)
        if gap > DX/2 + 1e-9:
            tail_ok = False
        #> under ord2 the exit row IS printed (x = L), so the closure term is ~0; it is kept
        #> as a guard in case the last row ever lands at the final dual face instead.
        dmdx = (rows[-1][8] - rows[-2][8]) / (rows[-1][0] - rows[-2][0]) if len(rows) > 1 else 0.0
        evap += npdot * ((m_first - m_last) + (-dmdx) * gap)

    resid = abs(src_total - evap) / evap if evap > 0 else 1.0
    n_exit = sum(1 for l in open(OUTLOC).read().splitlines()[2:] if len(l.split()) == 9)

    checks = [
        ("W1 gas-order = 2 was read (no ord1 fallback notice)",
         "Defaulting to 2nd order" not in log and "Output field: gas coupling source" in log),
        ("W2 mollify = off was read (witness)", "Field mollification OFF" in log),
        (f"P  all {NPART} parcels injected and exited", len(parts) == NPART and n_exit == NPART),
        (f"X  last row reaches the outlet: gap {worst_gap:.3e} m <= dx/2", tail_ok),
        (f"T  telescoping: sum(wdot)={src_total:.6e} vs evaporated={evap:.6e} kg/s, "
         f"resid={resid:.3e} <= {TELESCOPE_TOL}", evap > 0.0 and src_total > 0.0 and resid <= TELESCOPE_TOL),
    ]
    fail = 0
    for name, ok in checks:
        print(f"[{'PASS' if ok else 'FAIL'}] {name}")
        fail += (not ok)
    print(f"       mesh {nx}x{ny}x{nz}; parcels ride dual row j = 1 (a FACE row), "
          f"which the pre-O28 reduction handed on at 1/2 => resid ~ 0.50")
    return 1 if fail else 0


if __name__ == "__main__":
    sys.exit(main())
