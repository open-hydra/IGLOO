#!/usr/bin/env python3
"""Overlay IGLOO vs the Godsave-Spalding analytic d^2 line (E-VAL-1, MOSE-style).

check.py is the GATE (input-only closed-form line, theory budget); this only DRAWS
and writes OUTPUT/d2law-line.svg. matplotlib imported defensively (no-op if absent).

One panel: d^2/d0^2 vs x — the analytic line (gated reference, K0 from case inputs
alone) and IGLOO's recorded points; inset text reports K0 and the frozen Tp drift.
"""
import argparse
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
    ap = argparse.ArgumentParser(description="overlay IGLOO vs analytic d2 line")
    ap.add_argument("--plot", action="store_true")
    args = ap.parse_args()
    if not os.path.isfile(TRAJ):
        print(f"[verify] {TRAJ} not found -- run the case first")
        return 0
    parts = check.load(TRAJ)
    pid = next((p for p in sorted(parts) if len(parts[p]) >= check.MIN_PTS), None)
    if pid is None:
        print("[verify] no sufficiently sampled particle to plot")
        return 0
    rows = sorted(parts[pid])
    seen = [rows[0]]
    for r in rows[1:]:
        if r[0] > seen[-1][0]:
            seen.append(r)
    xa, _, da = seen[0]
    d0sq = da * da
    xs = [r[0] for r in seen]
    d2n = [(r[2] ** 2) / d0sq for r in seen]
    tp_drift = max(r[1] for r in seen) - check.TP0
    line = [(d0sq - check.K0 * (x - xa) / check.U_G) / d0sq for x in xs]

    fig, ax = plt.subplots(figsize=(6.2, 4.0))
    ax.plot(xs, line, "-", color="C0", lw=1.6,
            label="analytic line (gated): $d^2 = d_0^2 - K_0\\,x/u_g$")
    ax.plot(xs, d2n, "o", color="C1", ms=3.0, mfc="none", label="IGLOO")
    ax.set_xlabel("$x$ [m]"); ax.set_ylabel("$d^2/d_0^2$")
    ax.set_title(f"{CASE}: Godsave-Spalding $d^2$-law, input-only line "
                 f"[God53/Spa53]")
    ax.text(0.03, 0.08,
            f"$K_0$ = {check.K0:.4e} m$^2$/s (from inputs)\n"
            f"$B_{{T0}}$ = {check.BT0:.4f};  $T_p$ drift = {tp_drift:.2f} K "
            f"(cp$_p\\times$100 freeze)",
            transform=ax.transAxes, fontsize=7,
            bbox=dict(fc="white", alpha=0.7, ec="none"))
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
