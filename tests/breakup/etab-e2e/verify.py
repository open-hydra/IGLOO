#!/usr/bin/env python3
"""Overlay IGLOO vs Tanner ETAB (B-VAL-5, MOSE-style): child-size ratio vs We.

check.py is the GATE (onset + t_bu + cascade product size); this only DRAWS and writes
OUTPUT/etab-e2e.svg. matplotlib imported defensively.

One panel: d_child/d_parent vs We_r (log-x): the Tanner closed form (gated), the
bag/stripping branch split at WeTrans=80, the stripping asymptote exp(-k2*sqrt(24)),
and IGLOO's measured first-break ratios.
"""
import argparse
import math
import os
import sys

import check

HERE = os.path.dirname(os.path.abspath(__file__))
CASE = os.path.basename(HERE)
TRAJ = os.path.join(HERE, "OUTPUT", "trajectories-A.dat")

try:
    sys.path.insert(0, os.path.dirname(os.path.dirname(HERE)))
    import vv_style
    vv_style.apply()
    import matplotlib.pyplot as plt
except Exception as exc:
    print(f"[verify] plotting skipped ({exc.__class__.__name__}: {exc})")
    sys.exit(0)


def main():
    ap = argparse.ArgumentParser(description="overlay IGLOO vs Tanner ETAB size")
    ap.add_argument("--plot", action="store_true")
    args = ap.parse_args()
    if not os.path.isfile(TRAJ):
        print(f"[verify] {TRAJ} not found -- run the case first")
        return 0
    parts = check.load_trajectories(TRAJ)

    wes = [6.05 * (110.0 / 6.05) ** (i / 300.0) for i in range(301)]
    rts = [check.ratio_tanner(w) for w in wes]

    meas = []
    for pid in sorted(parts):
        raw = sorted(parts[pid])          # record order is OMP-nondeterministic
        rows = [raw[0]]
        for r in raw[1:]:
            if r[0] > rows[-1][0]:
                rows.append(r)
        if len(rows) < 4:
            continue
        _, tlo, _, _, we_brk, r_meas = check.drop_report(rows)
        if tlo is not None:
            meas.append((we_brk, r_meas))

    fig, ax = plt.subplots(figsize=(6.2, 4.0))
    ax.plot(wes, rts, "-", color="C0", lw=1.6, label="Tanner cascade (gated)")
    if meas:
        ax.plot(*zip(*meas), "o", color="C1", ms=4.5, label="IGLOO (measured 1st break)")
    asy = math.exp(-check.K2 * math.sqrt(24.0))
    ax.axhline(asy, ls=":", color="gray", lw=1,
               label=f"strip asymptote $e^{{-k_2\\sqrt{{24}}}}={asy:.3f}$")
    ax.axvline(check.WETRANS, ls="--", color="gray", lw=0.8, alpha=0.6)
    ax.text(check.WETRANS * 1.02, 0.72, "bag$\\,\\to\\,$strip", fontsize=7, color="gray")
    ax.set_xscale("log")
    ax.set_xlim(6, 110); ax.set_ylim(0.25, 0.85)
    ax.set_xlabel("$We_r = \\rho_g u^2 r/\\sigma$")
    ax.set_ylabel("$d_\\mathrm{child}/d_\\mathrm{parent}$")
    ax.set_title(f"{CASE}: ETAB cascade product size vs Tanner [Tan97/98]")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(loc="upper right")
    fig.tight_layout()
    out = os.path.join(HERE, "OUTPUT", f"{CASE}.svg")
    fig.savefig(out, bbox_inches="tight", transparent=True)
    print(f"[verify] wrote {out}")
    if args.plot:
        plt.show()
    return 0


if __name__ == "__main__":
    sys.exit(main())
