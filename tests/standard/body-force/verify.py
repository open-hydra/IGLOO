#!/usr/bin/env python3
"""Plot IGLOO body-force output vs the EXACT terminal-drift reference (MOSE-style).

check.py is the GATE; this only draws the overlay and writes OUTPUT/body-force.svg.
It never gates the suite. matplotlib is imported defensively (no-op if absent).

  ./verify.py [--plot]

Reference: reference/drift.txt (x, v), the exact transverse drift
  v(x) = v_inf (1 - exp(-x/L)),  v_inf = g_y tau,  tau = rho_p d^2 /(18 mu)
reproduced from the declared law in check.py (not digitized). The IGLOO transverse
velocity v(x) (col 5) is overlaid as markers.
"""
import argparse
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
CASE = os.path.basename(HERE)
TRAJ = os.path.join(HERE, "OUTPUT", "trajectories-A.dat")
REF  = os.path.join(HERE, "reference", "drift.txt")
VCOL, MIN_PTS = 4, 3                        # col 5 (0-based 4) = transverse velocity v

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
    ap = argparse.ArgumentParser(description="overlay IGLOO vs terminal-drift curve")
    ap.add_argument("--plot", action="store_true", help="also open a window")
    args = ap.parse_args()
    if not os.path.isfile(TRAJ):
        print(f"[verify] {TRAJ} not found -- run the case first")
        return 0
    ref = np.loadtxt(REF)
    xr, vr = ref[:, 0], ref[:, 1]
    got = best_particle(TRAJ, VCOL)
    if got is None:
        print("[verify] no sufficiently sampled particle to plot")
        return 0
    xa, xs, vs = got

    fig, ax = plt.subplots(figsize=(7.0, 4.2))
    ax.plot(xr + xa, vr, "-", color="C0", lw=1.6, label="body-force drift:  $v(x)=v_\\infty(1-e^{-x/L})$\n$v_\\infty=g_y\\tau$,  $L=u_g\\tau$,  $\\tau=\\rho_p d^2/(18\\mu)$", zorder=2)
    ax.plot(xs, vs, "o", color="C1", ms=3.5, mfc="none", label="IGLOO", zorder=3)
    ax.set_xlim(min(xs), max(xs))
    ax.set_xlabel("$x$ [m]")
    ax.set_ylabel("$v$ [m/s]")
    ax.set_title(f"{CASE}: uniform body-force terminal drift")
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
