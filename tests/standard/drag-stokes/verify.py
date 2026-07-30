#!/usr/bin/env python3
"""Plot IGLOO drag-stokes output vs the EXACT Stokes reference (MOSE-style verify.py).

Companion to check.py: check.py is the GATE (recomputes the oracle in-process and
asserts pass/fail); this script only draws the overlay and writes an SVG artifact.
It never gates the suite -- ctest runs check.py, not this file. matplotlib is
imported defensively so running this on a host without it prints a note and exits
0 instead of crashing.

  ./verify.py            # write OUTPUT/drag-stokes.svg
  ./verify.py --plot     # also open an interactive window

Reference curve: reference/stokes.txt (u, x), the exact inlet-anchored relaxation
  u(t) = u_g + (v0-u_g) e^{-t/tau},   x(u) = -u_g*tau*ln((u-u_g)/(v0-u_g)) - tau*(u-v0)
The IGLOO trajectory's axial velocity u(x) is overlaid as markers.
"""
import argparse
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
CASE = os.path.basename(HERE)
TRAJ = os.path.join(HERE, "OUTPUT", "trajectories-A.dat")
REF  = os.path.join(HERE, "reference", "stokes.txt")
V0, U_G = 1.0, 10.0
MIN_PTS = 3

try:
    import numpy as np
    sys.path.insert(0, os.path.dirname(os.path.dirname(HERE)))
    import vv_style                            # shared V&V style (CM/LaTeX math fonts)
    vv_style.apply()
    import matplotlib.pyplot as plt
except Exception as exc:                        # matplotlib/numpy absent -> no-op
    print(f"[verify] plotting skipped ({exc.__class__.__name__}: {exc})")
    sys.exit(0)


def load_trajectories(path):
    """{ID: [(x,u), ...]} axial position/velocity per particle, sorted by x."""
    parts = {}
    with open(path) as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith("variables") or s.startswith("Zone"):
                continue
            c = s.split()
            if len(c) < 10:
                continue
            try:
                x, u = float(c[0]), float(c[3])
                pid = int(c[9])
            except ValueError:
                continue
            parts.setdefault(pid, []).append((x, u))
    for pid in parts:
        parts[pid].sort(key=lambda r: r[0])
    return parts


def main():
    ap = argparse.ArgumentParser(description="overlay IGLOO vs exact Stokes curve")
    ap.add_argument("--plot", action="store_true", help="also open a window")
    args = ap.parse_args()

    if not os.path.isfile(TRAJ):
        print(f"[verify] {TRAJ} not found -- run the case first")
        return 0
    ref = np.loadtxt(REF)
    ur, xr = ref[:, 0], ref[:, 1]

    parts = load_trajectories(TRAJ)
    best = max((rows for rows in parts.values() if len(rows) >= MIN_PTS),
               key=len, default=None)
    if best is None:
        print("[verify] no sufficiently sampled particle to plot")
        return 0
    xa = best[0][0]
    xs = [r[0] for r in best]
    us = [r[1] for r in best]

    fig, ax = plt.subplots(figsize=(7.0, 4.2))
    ax.plot(xr + xa, ur, "-", color="C0", lw=1.6, label="exact Stokes:  $x(u)=x_a-u_g\\tau\\,\\ln\\frac{u-u_g}{u_a-u_g}-\\tau(u-u_a)$\n$\\tau=\\rho_p d^2/(18\\mu)$", zorder=2)
    ax.plot(xs, us, "o", color="C1", ms=3.5, mfc="none", label="IGLOO", zorder=3)
    ax.set_xlabel("$x$ [m]")
    ax.set_ylabel("$u$ [m/s]")
    ax.set_title(f"{CASE}: Stokes drag relaxation [SN33]")
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
