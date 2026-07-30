#!/usr/bin/env python3
"""Plot the db-injection DIAGNOSTIC (MOSE-style verify.py).

check.py is the GATE; this only draws the artifact and writes OUTPUT/db-injection.svg.
It never gates the suite. matplotlib is imported defensively (no-op if absent).

  ./verify.py [--plot]

No reference curve: this is an infrastructure case (assigned-position "DB"
injection, bug-B1 regression), not a physics-paper case. The plot shows the
axial velocity relaxation u(x) of the 5 assigned-position particles,
confirming the vInj hand-off (u starts at up=1) and Stokes relaxation.
"""
import argparse
import os
import sys

import check                                    # reuse the gate's loader

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


def main():
    ap = argparse.ArgumentParser(description="plot assigned-position u(x) diagnostic")
    ap.add_argument("--plot", action="store_true", help="also open a window")
    args = ap.parse_args()
    if not os.path.isfile(TRAJ):
        print(f"[verify] {TRAJ} not found -- run the case first")
        return 0
    traj = check.rows_of(TRAJ, 10)

    fig, ax = plt.subplots(figsize=(7.0, 4.2))
    for i, pid in enumerate(sorted(traj)):
        rows = traj[pid]
        xs, us = [rows[0][0]], [rows[0][3]]
        for r in rows[1:]:
            if r[0] > xs[-1]:
                xs.append(r[0])
                us.append(r[3])
        ax.plot(xs, us, "o-", color=f"C{i % 10}", ms=3.0, lw=1.0, mfc="none",
                label=f"ID {pid}")
    ax.set_xlabel("$x$ [m]")
    ax.set_ylabel("$u$ [m/s]")
    ax.set_title(f"{CASE}: assigned-position injection + Stokes relaxation (diagnostic)")
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
