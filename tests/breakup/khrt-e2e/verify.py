#!/usr/bin/env python3
"""Overlay IGLOO vs Reitz-87 KH stripping (B-VAL-6, MOSE-style).

check.py is the GATE (initial KH-stripping rate, We_r >= 340); this only DRAWS and writes
OUTPUT/khrt-e2e.svg. matplotlib imported defensively.

One panel: initial dd/dt vs We_r — the Reitz-87 KH-stripping rate (gated) and IGLOO's
measured initial rates. The RT/shed shatter (bug A19, fixed 2026-07-23) is gated
separately by check_rt.py; this panel is the continuous KH-stripping rate.
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
    ap = argparse.ArgumentParser(description="overlay IGLOO vs Reitz-87 KH stripping")
    ap.add_argument("--plot", action="store_true")
    args = ap.parse_args()
    if not os.path.isfile(TRAJ):
        print(f"[verify] {TRAJ} not found -- run the case first")
        return 0
    parts = check.load_trajectories(TRAJ)

    wes = [check.WELIM * 1.5 * (1400.0 / (check.WELIM * 1.5)) ** (i / 300.0)
           for i in range(301)]
    slip = 0.5 * check.U_GAS
    ref = []
    for we in wes:
        d = we * check.SIGMA / (check.RHOG * slip**2)   # d from We_r at this slip
        ref.append((we, check.kh_rate(d, slip)))
    ref = [(w, r) for w, r in ref if r < 0.0]

    gated, ungated = [], []
    for pid in sorted(parts):
        raw = sorted(parts[pid])          # record order is OMP-nondeterministic
        rows = [raw[0]]
        for r in raw[1:]:
            if r[0] > rows[-1][0]:
                rows.append(r)
        if len(rows) < 4:
            continue
        d0, u0 = rows[0][2], rows[0][1]
        wer = 0.5 * d0 * check.RHOG * (check.U_GAS - u0)**2 / check.SIGMA
        r_ig, _ = check.windowed_rates(rows)
        if r_ig is None:
            continue
        (gated if wer >= check.WE_MIN else ungated).append((wer, r_ig))

    fig, ax = plt.subplots(figsize=(6.2, 4.0))
    if ref:
        ax.plot(*zip(*ref), "-", color="C0", lw=1.6, label="Reitz-87 KH rate (gated)")
    if gated:
        ax.plot(*zip(*gated), "o", color="C1", ms=4.5, label="IGLOO (gated, RT-free window)")
    if ungated:
        ax.plot(*zip(*ungated), "o", color="C1", ms=4.5, mfc="none",
                label="IGLOO (RT would fire in-window; not gated)")
    ax.axvline(check.WE_MIN, ls=":", color="gray", lw=1)
    ax.set_xscale("log")
    ax.set_xlabel("$We_r = \\rho_g u^2 r/\\sigma$")
    ax.set_ylabel("initial $\\mathrm{d}d/\\mathrm{d}t$ [m/s]")
    ax.set_title(f"{CASE}: KH-stripping rate vs Reitz-87 [Reitz87] "
                 "(RT shatter gated by check_rt.py)")
    ax.grid(True, which="both", alpha=0.25)
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
