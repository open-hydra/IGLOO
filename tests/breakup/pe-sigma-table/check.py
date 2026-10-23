#!/usr/bin/env python3
"""Gate for pe-sigma-table: the breakup surface tension and viscosity read from properties.dat.

pilch-erdman-e2e's 25-drop Weber sweep (its gas file and boundary file) with sigma and mu taken from the
Sigma and Viscosity columns of INPUT/properties.dat, sigma = 0.036 + 1e-3 (T - 300) N/m and
mu = 0.05 exp(-(T - 300)/10) Pa s on 270..370 K, and no [IGLOO-Properties]. heat = NoHeat keeps every
drop at the injection temperature T0 = 300 K (the gas temperature), so the solver's sigma and mu are the
columns at 300 K: 0.036 N/m and 0.05 Pa s (half and fifty times the pilch-erdman-e2e constants).

Oracle: pilch-erdman-e2e's PE87 kernel (its check.py, imported unchanged) with its module constants SIGMA and
MUP set to the columns' values at T0, read from the file by tests/tools/proptab.py -- never through the
solver's reader. Gate (pilch-erdman-e2e's first gate): IGLOO's initial breakup rate of every drop with
We >= 45 within 1 % of the kernel's, the same-window least-squares slope on both sides.
"""
import importlib.util
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "tools"))
from proptab import load_properties, table_value  # noqa: E402

_spec = importlib.util.spec_from_file_location("pe_kernel", os.path.join(HERE, "..", "pilch-erdman-e2e", "check.py"))
pe = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(pe)

TRAJ = "OUTPUT/trajectories-A.dat"
T0 = 300.0
N_GOOD = 20       # gated drops: 24 on the reference run (drop 2, We 60, keeps its diameter on both sides)


def main():
    Tmin, Tmax, cols = load_properties("INPUT/properties.dat")
    pe.SIGMA = table_value(cols["sigma"], Tmin, Tmax, T0)
    pe.MUP = table_value(cols["mu"], Tmin, Tmax, T0)
    print(f"sigma({T0:.0f} K) = {pe.SIGMA:.6g} N/m, mu({T0:.0f} K) = {pe.MUP:.6g} Pa s (INPUT/properties.dat)")

    try:
        rows_T = [float(c.split()[6]) for c in open(TRAJ)
                  if len(c.split()) == 10 and c.split()[0][0] in "0123456789-"]
    except FileNotFoundError:
        print(f"[FAIL] {TRAJ} not found -- did the solver run?")
        return 1
    off = [t for t in rows_T if t != T0]
    print(f"trajectory rows: {len(rows_T)}, rows off T0: {len(off)}")
    if not rows_T or off:
        print("[FAIL] the drop temperature is not T0 on every row: the oracle's sigma and mu do not apply")
        return 1

    parts = pe.load_trajectories(TRAJ)
    print(f"{'ID':>3} {'We0':>6} {'rate_IG':>10} {'rate_PE87':>10} {'ratio':>6}  verdict")
    n_good = n_viol = 0
    worst = 0.0
    for pid in sorted(parts):
        raw = sorted(parts[pid])          # record order is OMP-nondeterministic
        rows = [raw[0]]
        for r in raw[1:]:
            if r[0] > rows[-1][0]:
                rows.append(r)
        if len(rows) < 4:
            continue
        d0, u0 = rows[0][2], rows[0][1]
        slip0 = pe.U_GAS - u0
        we0 = d0 * pe.RHOG * slip0**2 / pe.SIGMA
        if we0 < pe.WE_MIN_GATE:
            continue
        r_ig, r_pe = pe.windowed_rates(rows)
        if r_ig is None or r_pe is None or r_pe >= 0.0:
            continue
        rel = abs(r_ig - r_pe) / abs(r_pe)
        n_good += 1
        worst = max(worst, rel)
        ok = rel <= pe.TOL_REL
        if not ok:
            n_viol += 1
        print(f"{pid:>3} {we0:>6.0f} {r_ig:>10.3e} {r_pe:>10.3e} {r_ig / r_pe:>6.3f}  {'PASS' if ok else 'FAIL'}")

    print(f"\ndrops gated (We >= {pe.WE_MIN_GATE:.0f}): {n_good} (need >= {N_GOOD}); "
          f"violations (rel > {pe.TOL_REL:.0%}): {n_viol}; worst rel {worst:.2e}")
    if n_good >= N_GOOD and n_viol == 0:
        print("\n[PASS] the PE87 rate at the table's sigma and mu at every gated drop (tol 1%).")
        return 0
    print(f"\n[FAIL] gated drops {n_good} (need >= {N_GOOD}), rate violations {n_viol}.")
    return 1


if __name__ == "__main__":
    sys.exit(main())
