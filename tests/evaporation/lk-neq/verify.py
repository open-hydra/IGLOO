#!/usr/bin/env python3
"""Plot IGLOO lk-neq output vs the Langmuir-Knudsen law (MOSE-style verify.py).

check.py is the GATE; this only draws the overlay and writes OUTPUT/lk-neq.svg.
It never gates the suite. matplotlib is imported defensively (no-op if absent).

  ./verify.py [--plot]

Reference: the declared LK non-equilibrium law [MHB98] is REPRODUCED (not
digitized) from the gate's own oracle in check.py -- imported here as a single
source of truth. The rate law is an ODE in d² driven by the particle temperature,
so the reference d(x) is integrated along IGLOO's MEASURED Tp (col 7), exactly as
the gate does; there is no static reference/ data file (it is run-conditioned).
The VLE equilibrium curve is drawn dashed as the regression the case out-distances.
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
    """(pid, xs, d_meas_um, d_lk_um, d_vle_um) of the strongest evaporator."""
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
            pv.append(check.rk4_span(xs[k-1], xs[k], Tps[k-1], Tps[k], pv[-1], equilibrium=True))
        rel_loss = (d2a - rows[-1][7] ** 2) / d2a
        cand = (rel_loss, pid, xs,
                [r[7] * 1e6 for r in rows],
                [math.sqrt(max(v, 0.0)) * 1e6 for v in pf],
                [math.sqrt(max(v, 0.0)) * 1e6 for v in pv])
        if best is None or rel_loss > best[0]:
            best = cand
    return None if best is None else best[1:]


def main():
    ap = argparse.ArgumentParser(description="overlay IGLOO vs Langmuir-Knudsen law")
    ap.add_argument("--plot", action="store_true", help="also open a window")
    args = ap.parse_args()
    if not os.path.isfile(TRAJ):
        print(f"[verify] {TRAJ} not found -- run the case first")
        return 0
    got = best_evaporator(check.load_trajectories(TRAJ))
    if got is None:
        print("[verify] no sufficiently sampled evaporator to plot")
        return 0
    pid, xs, d_meas, d_lk, d_vle = got

    fig, ax = plt.subplots(figsize=(7.0, 4.2))
    ax.plot(xs, d_lk, "-", color="C0", lw=1.6, label="LK non-eq.:  $\\frac{d(d^2)}{dt}=-\\frac{8\\rho_g D_v}{\\rho_p}\\ln(1+B_M)$\n$X_s=X_s^{eq}-\\frac{2L_K}{d}\\beta$,  $\\beta=\\frac{-\\dot{m}\\,\\mathrm{Pr}}{2\\pi\\mu_g d}$", zorder=2)
    ax.plot(xs, d_vle, "--", color="C3", lw=1.2, label="VLE equilibrium (regression)", zorder=1)
    ax.plot(xs, d_meas, "o", color="C1", ms=3.5, mfc="none", label="IGLOO", zorder=3)
    ax.set_xlabel("$x$ [m]")
    ax.set_ylabel("$d$ [$\\mu$m]")
    ax.set_title(f"{CASE}: Langmuir-Knudsen non-equilibrium along measured $T_p$ [MHB98] (ID {pid})")
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
