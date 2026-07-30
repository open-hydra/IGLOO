#!/usr/bin/env python3
"""Plot the coupled-body transverse drift vs the closed form (MOSE-style overlay).

check.py is the GATE; this only draws the v(x) overlay (same physics as
standard/body-force, with euler+source+body all active) and writes
OUTPUT/coupled-body.svg. matplotlib is imported defensively (no-op if absent).

Reference (inline, not digitized): v(x) = v_inf (1 - e^(-x/L)),
v_inf = g_y*tau, L = u_g*tau, tau = rho_p d^2/(18 mu).
"""
import argparse
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
CASE = os.path.basename(HERE)
TRAJ = os.path.join(HERE, "OUTPUT", "trajectories-A.dat")

RHO_P, D_P, MU = 2950.0, 1.189e-5, 1.8e-5
U_G, G_Y = 10.0, -200.0
TAU = RHO_P * D_P**2 / (18.0 * MU)
L, V_INF = U_G * TAU, G_Y * TAU
VCOL, MIN_PTS = 4, 3

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
            c = line.split()
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
    ap.add_argument("--plot", action="store_true")
    args = ap.parse_args()
    if not os.path.isfile(TRAJ):
        print(f"[verify] {TRAJ} not found -- run the case first")
        return 0
    got = best_particle(TRAJ, VCOL)
    if got is None:
        print("[verify] no sufficiently sampled particle to plot")
        return 0
    xa, xs, vs = got
    xg = np.linspace(min(xs), max(xs), 400)
    vg = V_INF * (1.0 - np.exp(-(xg - xa) / L))

    fig, ax = plt.subplots(figsize=(7.0, 4.2))
    ax.plot(xg, vg, "-", color="C0", lw=1.6,
            label="body-force drift (euler+source+body on):\n$v(x)=v_\\infty(1-e^{-x/L})$,  $v_\\infty=g_y\\tau$,  $L=u_g\\tau$", zorder=2)
    ax.plot(xs, vs, "o", color="C1", ms=3.5, mfc="none", label="IGLOO", zorder=3)
    ax.set_xlabel("$x$ [m]")
    ax.set_ylabel("$v$ [m/s]")
    ax.set_title(f"{CASE}: combined accumulators, trajectory untouched")
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
