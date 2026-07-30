#!/usr/bin/env python3
"""Plot the db-2daxi DIAGNOSTIC (MOSE-style verify.py).

check.py is the GATE; this only draws the artifact and writes OUTPUT/db-2daxi.svg.
It never gates the suite. matplotlib is imported defensively (no-op if absent).

  ./verify.py [--plot]

No reference curve: this is an infrastructure case (axisymmetric-wedge mesh,
real MOSE nozzle flow, DB injection, euler output on), not a physics-paper
case. Two panels track both particles through the nozzle: the meridional path
y(x) (delthe fold keeps particles in the wedge) and the temperature history
T_p(x) through the ~300-3600 K gas field, throat at x = 0.175 m.
"""
import argparse
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
CASE = os.path.basename(HERE)
TRAJ = os.path.join(HERE, "OUTPUT", "trajectories-A.dat")
X_THROAT = 0.175

try:
    sys.path.insert(0, os.path.dirname(os.path.dirname(HERE)))
    import vv_style                            # shared V&V style (CM/LaTeX math fonts)
    vv_style.apply()
    import matplotlib.pyplot as plt
except Exception as exc:
    print(f"[verify] plotting skipped ({exc.__class__.__name__}: {exc})")
    sys.exit(0)


def load_trajectories(path):
    """{ID: [(x, y, T), ...]} in file order (rows are already time-ordered)."""
    parts = {}
    with open(path) as f:
        for line in f:
            t = line.split()
            if len(t) != 10:
                continue
            try:
                pid = int(t[-1])
                x, y, T = float(t[0]), float(t[1]), float(t[6])
            except ValueError:
                continue
            parts.setdefault(pid, []).append((x, y, T))
    return parts


def main():
    ap = argparse.ArgumentParser(description="plot db-2daxi nozzle-transit diagnostic")
    ap.add_argument("--plot", action="store_true", help="also open a window")
    args = ap.parse_args()
    if not os.path.isfile(TRAJ):
        print(f"[verify] {TRAJ} not found -- run the case first")
        return 0
    parts = load_trajectories(TRAJ)
    if not parts:
        print("[verify] no trajectory rows to plot")
        return 0

    fig, axs = plt.subplots(2, 1, figsize=(7.0, 6.4), sharex=True)
    for i, pid in enumerate(sorted(parts)):
        rows = parts[pid]
        xs = [r[0] for r in rows]
        ys = [r[1] for r in rows]
        Ts = [r[2] for r in rows]
        axs[0].plot(xs, ys, "o-", color=f"C{i % 10}", ms=2.5, lw=0.9, mfc="none",
                    label=f"ID {pid}")
        axs[1].plot(xs, Ts, "o-", color=f"C{i % 10}", ms=2.5, lw=0.9, mfc="none",
                    label=f"ID {pid}")
    for ax in axs:
        ax.axvline(X_THROAT, color="0.5", lw=0.8, ls=":")
        ax.grid(True, alpha=0.25)
        ax.legend(loc="best")
    axs[0].set_ylabel("$y$ [m]")
    axs[0].set_title(f"{CASE}: 2D-axisymmetric nozzle transit (diagnostic)")
    axs[1].set_ylabel("$T_p$ [K]")
    axs[1].set_title("particle temperature through the nozzle (dotted: throat)")
    axs[1].set_xlabel("$x$ [m]")
    fig.tight_layout()
    out = os.path.join(HERE, "OUTPUT", f"{CASE}.svg")
    fig.savefig(out, bbox_inches="tight", transparent=True)
    print(f"[verify] wrote {out}")
    if args.plot:
        plt.show()
    return 0


if __name__ == "__main__":
    sys.exit(main())
