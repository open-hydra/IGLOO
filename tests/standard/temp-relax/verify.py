#!/usr/bin/env python3
"""Plot IGLOO temp-relax output vs the EXACT Nu=2 reference (MOSE-style verify.py).

check.py is the GATE; this only draws the overlay and writes OUTPUT/temp-relax.svg.
It never gates the suite. matplotlib is imported defensively (no-op if absent).

  ./verify.py [--plot]

Reference: reference/nu2.txt (x, T), the exact lumped-capacitance relaxation
  T(x) = T_g + (Tp0-T_g) exp(-x/L),  L = u_g cp rho_p d^2 /(12 k_g)   (Nu=2)
reproduced from the declared law in check.py (not digitized). The IGLOO
trajectory temperature T_p(x) (col 7) is overlaid as markers.
"""
import argparse
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
CASE = os.path.basename(HERE)
TRAJ = os.path.join(HERE, "OUTPUT", "trajectories-A.dat")
REF  = os.path.join(HERE, "reference", "nu2.txt")
VCOL, MIN_PTS = 6, 3                       # col 7 (0-based 6) = temperature

try:
    import numpy as np
    sys.path.insert(0, os.path.dirname(os.path.dirname(HERE)))
    import vv_style                            # shared V&V style (CM/LaTeX math fonts)
    vv_style.apply()
    import matplotlib.pyplot as plt
except Exception as exc:
    print(f"[verify] plotting skipped ({exc.__class__.__name__}: {exc})")
    sys.exit(0)


def best_particle(path, vcol):
    """(xa, xs, ys) of the best-sampled particle: x and column `vcol`."""
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
                x, y = float(c[0]), float(c[vcol])
                pid = int(c[9])
            except ValueError:
                continue
            parts.setdefault(pid, []).append((x, y))
    best = max((sorted(r) for r in parts.values() if len(r) >= MIN_PTS),
               key=len, default=None)
    if best is None:
        return None
    return best[0][0], [r[0] for r in best], [r[1] for r in best]


def main():
    ap = argparse.ArgumentParser(description="overlay IGLOO vs exact Nu=2 curve")
    ap.add_argument("--plot", action="store_true", help="also open a window")
    args = ap.parse_args()
    if not os.path.isfile(TRAJ):
        print(f"[verify] {TRAJ} not found -- run the case first")
        return 0
    ref = np.loadtxt(REF)
    xr, Tr = ref[:, 0], ref[:, 1]
    got = best_particle(TRAJ, VCOL)
    if got is None:
        print("[verify] no sufficiently sampled particle to plot")
        return 0
    xa, xs, Ts = got

    fig, ax = plt.subplots(figsize=(7.0, 4.2))
    ax.plot(xr + xa, Tr, "-", color="C0", lw=1.6, label="exact $\\mathrm{Nu}=2$:  $T_p(x)=T_g+(T_{p0}-T_g)\\,e^{-x/L}$\n$L=u_g\\rho_p c_p d^2/(12k_g)$", zorder=2)
    ax.plot(xs, Ts, "o", color="C1", ms=3.5, mfc="none", label="IGLOO", zorder=3)
    ax.set_xlim(min(xs), max(xs))
    ax.set_xlabel("$x$ [m]")
    ax.set_ylabel("$T_p$ [K]")
    ax.set_title(f"{CASE}: $\\mathrm{{Nu}}=2$ lumped-capacitance relaxation [RM52]")
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
