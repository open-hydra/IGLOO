#!/usr/bin/env python3
"""Plot IGLOO solid-box output vs the piecewise closed form (MOSE-style verify.py).

check.py is the GATE; this only draws the overlay and writes OUTPUT/solid-box.svg.
It never gates the suite. matplotlib is imported defensively (no-op if absent).

  ./verify.py [--plot]

Reference: the closed form of check.py -- liquid relaxation to the nucleation point x_n, the
recalescence to T_m, the freezing plateau at T_m to x_s, then solid relaxation with c_s.
"""
import argparse
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
CASE = os.path.basename(HERE)
sys.path.insert(0, HERE)
import check                                    # reuse the gate's constants and loader

try:
    import numpy as np
    sys.path.insert(0, os.path.dirname(os.path.dirname(HERE)))
    import vv_style                            # shared V&V style (CM/LaTeX math fonts)
    vv_style.apply()
    import matplotlib.pyplot as plt
except Exception as exc:
    print(f"[verify] plotting skipped ({exc.__class__.__name__}: {exc})")
    sys.exit(0)


def closed_form(x):
    if x < check.X_N:
        return check.T_G + (check.T_0 - check.T_G) * math.exp(-x / check.L_L)
    if x < check.X_S:
        return check.T_M
    return check.T_G + (check.T_M - check.T_G) * math.exp(-(x - check.X_S) / check.L_S)


def main():
    ap = argparse.ArgumentParser(description="overlay IGLOO vs the solidification closed form")
    ap.add_argument("--plot", action="store_true", help="also open a window")
    args = ap.parse_args()
    traj = os.path.join(HERE, check.TRAJ)
    if not os.path.isfile(traj):
        print(f"[verify] {traj} not found -- run the case first")
        return 0
    parts = check.load_rows(traj, 10)
    if not parts:
        print("[verify] no trajectory rows to plot")
        return 0
    rows = sorted(max(parts.values(), key=len), key=lambda r: r[0])
    xs, Ts = [r[0] for r in rows], [r[6] for r in rows]
    xr = np.linspace(0.0, check.LX, 1201)
    Tr = [closed_form(x) for x in xr]

    fig, ax = plt.subplots(figsize=(7.0, 4.2))
    ax.plot(xr, Tr, "-", color="C0", lw=1.6, zorder=2,
            label="closed form: liquid $L_l$, recalescence at $x_n$,\nplateau at $T_m$ to $x_s$, solid $L_s$")
    ax.plot(xs, Ts, "o", color="C1", ms=3.5, mfc="none", label="IGLOO", zorder=3)
    ax.axhline(check.T_N, color="0.6", lw=0.8, ls=":", zorder=1)
    ax.set_xlim(0.0, check.LX)
    ax.set_xlabel("$x$ [m]")
    ax.set_ylabel("$T_p$ [K]")
    ax.set_title(f"{CASE}: supercooling, recalescence and freezing plateau")
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
