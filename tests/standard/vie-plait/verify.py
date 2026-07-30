#!/usr/bin/env python3
"""Plot the Vié "plait": IGLOO trajectories vs the analytic damped oscillator.

check.py is the GATE; this only draws the overlay and writes OUTPUT/vie-plait.svg.
It never gates the suite. matplotlib is imported defensively (no-op if absent).

  ./verify.py [--plot]

Each particle injected into the compressive field v_g=-eps*(y-1) traces a decaying
transverse oscillation (Eq. 5.3); the strands cross at the axis y=1 (particle
trajectory crossing) forming the braid/"plait". Markers = IGLOO; solid lines =
the anchored closed form the gate checks against (reuses check.py's oracle).
"""
import argparse
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
CASE = os.path.basename(HERE)
TRAJ = os.path.join(HERE, "OUTPUT", "trajectories-A.dat")

try:
    sys.path.insert(0, os.path.dirname(os.path.dirname(HERE)))
    import numpy as np
    import vv_style                            # shared V&V style (CM/LaTeX math)
    vv_style.apply()
    import matplotlib.pyplot as plt
    import check                               # reuse the gate oracle + constants
except Exception as exc:
    print(f"[verify] plotting skipped ({exc.__class__.__name__}: {exc})")
    sys.exit(0)


def main():
    ap = argparse.ArgumentParser(description="plot the Vie compressive-field plait")
    ap.add_argument("--plot", action="store_true", help="also open a window")
    args = ap.parse_args()
    if not os.path.isfile(TRAJ):
        print(f"[verify] {TRAJ} not found -- run the case first")
        return 0
    parts = check.load(TRAJ)

    fig, ax = plt.subplots(figsize=(7.6, 4.2))
    ax.axhline(check.Y_AXIS, color="0.6", lw=0.7, ls=":")
    for i, pid in enumerate(sorted(parts)):
        rows = parts[pid]
        xa, ya, va = rows[0][0], rows[0][1], rows[0][2]
        Ya = ya - check.Y_AXIS
        xs = [r[0] for r in rows]
        ys = [r[1] for r in rows]
        yo = [check.y_oracle((x - xa) / check.UG0, Ya, va) for x in xs]
        c = f"C{i % 10}"
        ax.plot(xs, yo, "-", color=c, lw=1.3, zorder=2)
        ax.plot(xs, ys, "o", color=c, ms=2.2, mfc="none", zorder=3,
                label=f"$y_0={ya:.1f}$")
    ax.set_xlabel("$x$ [m]")
    ax.set_ylabel("$y$ [m]")
    ax.set_title(rf"{CASE}: inertial particles in $v_g=-\epsilon(y-1)$, "
                 rf"$\mathrm{{St}}={check.EPS*check.TAU:.0f}$ [Vie15]")
    ax.set_xlim(0, check.LX)
    ax.grid(True, alpha=0.25)
    # one proxy entry explaining marker vs line, plus the strands
    ax.plot([], [], "-", color="0.3", label="analytic Eq. 5.4")
    ax.plot([], [], "o", color="0.3", mfc="none", label="IGLOO")
    ax.legend(loc="upper right", ncol=2)
    fig.tight_layout()
    out = os.path.join(HERE, "OUTPUT", f"{CASE}.svg")
    fig.savefig(out, bbox_inches="tight", transparent=True)
    print(f"[verify] wrote {out}")
    if args.plot:
        plt.show()
    return 0


if __name__ == "__main__":
    sys.exit(main())
