#!/usr/bin/env python3
"""Plot IGLOO tc-box output vs the Tonini-Cossali law (MOSE-style verify.py).

check.py is the GATE; this only draws the overlay and writes OUTPUT/tc-box.svg.
It never gates the suite. matplotlib is imported defensively (no-op if absent).

  ./verify.py [--plot]

Reference: the declared Tonini-Cossali Stefan-Fuchs law [TC2012] is REPRODUCED
(not digitized) from the gate's own oracle in check.py -- imported here as a single
source of truth. The rate law is an ODE in d² driven by the particle temperature,
so the reference d(x) is integrated along IGLOO's MEASURED Tp (col 7), exactly as
the gate does; there is no static reference/ data file (it is run-conditioned). The
CEM Spalding curve is drawn dashed as the silent regression the case out-distances.
"""
import argparse
import math
import os
import sys

import check                                    # reuse the gate's oracle

HERE = os.path.dirname(os.path.abspath(__file__))
CASE = os.path.basename(HERE)
TRAJ = os.path.join(HERE, "OUTPUT", "trajectories-A.dat")

try:
    sys.path.insert(0, os.path.dirname(os.path.dirname(HERE)))
    import vv_style                            # shared V&V style (CM/LaTeX math fonts)
    vv_style.apply()
    import matplotlib.pyplot as plt
except Exception as exc:
    print(f"[verify] plotting skipped ({exc.__class__.__name__}: {exc})")
    sys.exit(0)


def best_evaporator(parts):
    """(pid, xs, d_meas_um, d_tc_um, d_cem_um) of the strongest evaporator."""
    best = None
    for pid in sorted(parts):
        raw = sorted(parts[pid])          # record order is OMP-nondeterministic
        rows = [raw[0]]
        for r in raw[1:]:
            if r[0] > rows[-1][0]:
                rows.append(r)
        if len(rows) <= check.MIN_PTS:
            continue
        xs = [r[0] for r in rows]
        Tps = [r[6] for r in rows]
        d2a = rows[0][7] ** 2
        pf, pv = [d2a], [d2a]
        for k in range(1, len(rows)):
            pf.append(check.rk4_span(xs[k-1], xs[k], Tps[k-1], Tps[k], pf[-1]))
            pv.append(check.rk4_span(xs[k-1], xs[k], Tps[k-1], Tps[k], pv[-1], cem=True))
        rel_loss = (d2a - rows[-1][7] ** 2) / d2a
        cand = (rel_loss, pid, xs,
                [r[7] * 1e6 for r in rows],
                [math.sqrt(max(v, 0.0)) * 1e6 for v in pf],
                [math.sqrt(max(v, 0.0)) * 1e6 for v in pv])
        if best is None or rel_loss > best[0]:
            best = cand
    return None if best is None else best[1:]


def main():
    ap = argparse.ArgumentParser(description="overlay IGLOO vs Tonini-Cossali law")
    ap.add_argument("--plot", action="store_true", help="also open a window")
    args = ap.parse_args()
    if not os.path.isfile(TRAJ):
        print(f"[verify] {TRAJ} not found -- run the case first")
        return 0
    got = best_evaporator(check.load_trajectories(TRAJ))
    if got is None:
        print("[verify] no sufficiently sampled evaporator to plot")
        return 0
    pid, xs, d_meas, d_tc, d_cem = got

    fig, ax = plt.subplots(figsize=(7.0, 4.2))
    ax.plot(xs, d_tc, "-", color="C0", lw=1.6, label="Tonini-Cossali:  $\\frac{d(d^2)}{dt}=-\\frac{8\\rho_g D_v}{\\rho_p}\\,\\hat{m}$\n$\\hat{m}+(T_s/T_g-1)\\,Le_v[f(\\hat{m}/Le_v)-1]=\\frac{M_v}{M_g}\\ln\\frac{1-X_\\infty}{1-X_s}$\n$f(x)=x/(1-e^{-x})$", zorder=2)
    ax.plot(xs, d_cem, "--", color="C3", lw=1.2, label="CEM Spalding (regression)", zorder=1)
    ax.plot(xs, d_meas, "o", color="C1", ms=3.5, mfc="none", label="IGLOO", zorder=3)
    ax.set_xlabel("$x$ [m]")
    ax.set_ylabel("$d$ [$\\mu$m]")
    ax.set_title(f"{CASE}: Tonini-Cossali analytical evaporation along measured $T_p$ [TC2012] (ID {pid})")
    ax.grid(True, alpha=0.25)
    ax.legend(loc="best")
    fig.tight_layout()
    out = os.path.join(HERE, "OUTPUT", f"{CASE}.svg")
    fig.savefig(out, bbox_inches="tight", transparent=True)
    print(f"[verify] wrote {out}")
    if args.plot:
        plt.show()
    return 0


if __name__ == "__main__":
    sys.exit(main())
