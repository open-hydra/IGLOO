#!/usr/bin/env python3
"""Plot IGLOO d2law output vs the Godsave-Spalding d²-law (MOSE-style verify.py).

check.py is the GATE; this only draws the overlay and writes OUTPUT/d2law.svg.
It never gates the suite. matplotlib is imported defensively (no-op if absent).

  ./verify.py [--plot]

Reference: the declared d²-law is REPRODUCED (not digitized) from the gate's own
oracle in check.py -- imported here as a single source of truth. Because the rate
kernel depends on the particle temperature, the reference d(x) is the kernel
integrated along IGLOO's MEASURED Tp (col 7), exactly as the gate does; there is
therefore no static reference/ data file (the reference is run-conditioned). The
strongest evaporator is shown, d (col 8) overlaid as markers vs the kernel.
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
    """(pid, xs, d_meas_um, d_pred_um) of the strongest evaporator."""
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
        d2a = rows[0][7] ** 2
        g = [check.integrand(r[6]) for r in rows]
        I, d_meas, d_pred = 0.0, [rows[0][7] * 1e6], [rows[0][7] * 1e6]
        for k in range(1, len(rows)):
            I += 0.5 * (g[k - 1] + g[k]) * (xs[k] - xs[k - 1]) / check.U_G
            d2p = d2a - I
            d_pred.append(math.sqrt(d2p) * 1e6 if d2p > 0.0 else 0.0)
            d_meas.append(rows[k][7] * 1e6)
        rel_loss = (d2a - rows[-1][7] ** 2) / d2a
        if best is None or rel_loss > best[0]:
            best = (rel_loss, pid, xs, d_meas, d_pred)
    return None if best is None else best[1:]


def main():
    ap = argparse.ArgumentParser(description="overlay IGLOO vs Godsave-Spalding kernel")
    ap.add_argument("--plot", action="store_true", help="also open a window")
    args = ap.parse_args()
    if not os.path.isfile(TRAJ):
        print(f"[verify] {TRAJ} not found -- run the case first")
        return 0
    got = best_evaporator(check.load_trajectories(TRAJ))
    if got is None:
        print("[verify] no sufficiently sampled evaporator to plot")
        return 0
    pid, xs, d_meas, d_pred = got

    fig, ax = plt.subplots(figsize=(7.0, 4.2))
    ax.plot(xs, d_pred, "-", color="C0", lw=1.6, label="Godsave-Spalding:  $\\frac{d(d^2)}{dt}=-\\frac{8k_g}{c_{p,g}\\rho_p}\\ln(1+B_T)$\n$B_T=c_{p,g}(T_g-T_p)/L_v$", zorder=2)
    ax.plot(xs, d_meas, "o", color="C1", ms=3.5, mfc="none", label="IGLOO", zorder=3)
    ax.set_xlabel("$x$ [m]")
    ax.set_ylabel("$d$ [$\\mu$m]")
    ax.set_title(f"{CASE}: Godsave-Spalding $d^2$-law along measured $T_p$ [God53/Spa53] (ID {pid})")
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
