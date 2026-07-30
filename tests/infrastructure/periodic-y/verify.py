#!/usr/bin/env python3
"""Plot the wrapped periodic-y trajectory vs the closed form (MOSE-style overlay).

check.py is the GATE; this only draws y(x) for the best-sampled particle with the
unwrapped body-force closed form folded into [0, Ly), and writes
OUTPUT/periodic-y.svg. matplotlib is imported defensively (no-op if absent).

Reference (inline, not digitized): y(x) = ya + v_inf*dt + (va - v_inf)*tau*(1-e^(-dt/tau)),
dt = (x-xa)/u_g, folded modulo Ly — the 201 transport is an exact +-Ly translation.
"""
import argparse
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
CASE = os.path.basename(HERE)
TRAJ = os.path.join(HERE, "OUTPUT", "trajectories-A.dat")

RHO_P, D_P, MU = 2950.0, 1.189e-5, 1.8e-5
U_G, G_Y, LY = 10.0, -8000.0, 0.05
TAU = RHO_P * D_P**2 / (18.0 * MU)
V_INF = G_Y * TAU
MIN_PTS = 10

try:
    import numpy as np
    sys.path.insert(0, os.path.dirname(os.path.dirname(HERE)))
    import vv_style                            # shared V&V style (CM/LaTeX math fonts)
    vv_style.apply()
    import matplotlib.pyplot as plt
except Exception as exc:
    print(f"[verify] plotting skipped ({exc.__class__.__name__}: {exc})")
    sys.exit(0)


def best_particle(path):
    parts = {}
    with open(path) as f:
        for line in f:
            c = line.split()
            if len(c) < 10:
                continue
            try:
                x, y, v = float(c[0]), float(c[1]), float(c[4])
                pid = int(c[9])
            except ValueError:
                continue
            parts.setdefault(pid, []).append((x, y, v))
    best = max((sorted(r) for r in parts.values() if len(r) >= MIN_PTS),
               key=len, default=None)
    return best


def main():
    ap = argparse.ArgumentParser(description="overlay wrapped y(x) vs closed form")
    ap.add_argument("--plot", action="store_true")
    args = ap.parse_args()
    if not os.path.isfile(TRAJ):
        print(f"[verify] {TRAJ} not found -- run the case first")
        return 0
    rows = best_particle(TRAJ)
    if rows is None:
        print("[verify] no sufficiently sampled particle to plot")
        return 0
    xa, ya, va = rows[0]
    xs = np.array([r[0] for r in rows])
    ys = np.array([r[1] for r in rows])

    xg = np.linspace(xs.min(), xs.max(), 800)
    dt = (xg - xa) / U_G
    y_unwrapped = ya + V_INF * dt + (va - V_INF) * TAU * (1.0 - np.exp(-dt / TAU))
    y_folded = np.mod(y_unwrapped, LY)
    # break the reference line at the wraps so matplotlib doesn't draw verticals
    jump = np.abs(np.diff(y_folded)) > 0.5 * LY
    y_plot = y_folded.copy()
    y_plot[1:][jump] = np.nan

    fig, ax = plt.subplots(figsize=(7.0, 4.2))
    ax.plot(xg, y_plot, "-", color="C0", lw=1.6,
            label="closed form mod $L_y$ (201 = exact $\\pm L_y$):\n$y=y_a+v_\\infty\\Delta t+(v_a-v_\\infty)\\,\\tau(1-e^{-\\Delta t/\\tau})$,  $\\Delta t=(x-x_a)/u_g$", zorder=2)
    ax.plot(xs, ys, "o", color="C1", ms=3.5, mfc="none", label="IGLOO", zorder=3)
    ax.axhline(0.0, color="0.5", lw=0.8, ls=":")
    ax.axhline(LY, color="0.5", lw=0.8, ls=":")
    ax.set_xlabel("$x$ [m]")
    ax.set_ylabel("$y$ [m]")
    ax.set_title(f"{CASE}: translational periodic transport (bcdef 201)")
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
