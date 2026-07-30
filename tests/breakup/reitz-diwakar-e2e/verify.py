#!/usr/bin/env python3
"""Overlay IGLOO vs Reitz-Diwakar [RD87] breakup rate (B-VAL-4, MOSE-style).

check.py is the GATE (initial breakup rate, bag+stripping); this only DRAWS and writes
OUTPUT/reitz-diwakar-e2e.svg. matplotlib imported defensively.

One panel: initial dd/dt vs We_r — the RD closed-form rate (gated), IGLOO's measured
initial rates coloured by regime (bag vs stripping), and the bag→stripping handoff line.
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
    ap = argparse.ArgumentParser(description="overlay IGLOO vs RD rate")
    ap.add_argument("--plot", action="store_true")
    args = ap.parse_args()
    if not os.path.isfile(TRAJ):
        print(f"[verify] {TRAJ} not found -- run the case first")
        return 0
    parts = check.load_trajectories(TRAJ)

    # analytic RD initial rate at injection (d0, slip0=U/2 by design) over We_r
    slip = 0.5 * check.U_GAS
    wes, rr = [], []
    for i in range(400):
        wer = 6.5 * (1050.0 / 6.5) ** (i / 399.0)
        d = 2.0 * wer * check.SIGMA / (check.RHOG * slip**2)   # We_r=ρg u² r/σ -> d
        wes.append(wer); rr.append(check.rd_rate(d, slip))

    bag, strip = [], []
    for pid in sorted(parts):
        raw = sorted(parts[pid])          # record order is OMP-nondeterministic
        rows = [raw[0]]
        for r in raw[1:]:
            if r[0] > rows[-1][0]:
                rows.append(r)
        if len(rows) < 4:
            continue
        d0, u0 = rows[0][2], rows[0][1]
        slip0 = check.U_GAS - u0
        wer = 0.5 * d0 * check.RHOG * slip0**2 / check.SIGMA
        r_ig, _ = check.windowed_rates(rows)
        if r_ig is None:
            continue
        (bag if check.rd_regime(d0, slip0) == "bag" else strip).append((wer, r_ig))

    fig, ax = plt.subplots(figsize=(6.2, 4.0))
    ax.plot(wes, rr, "-", color="C0", lw=1.6, label="RD [RD87] rate (gated)")
    if bag:
        ax.plot(*zip(*bag), "o", color="C1", ms=4.5, label="IGLOO (bag)")
    if strip:
        ax.plot(*zip(*strip), "s", color="C2", ms=4.0, label="IGLOO (stripping)")
    ax.axvline(20.0, ls="--", color="gray", lw=0.8, alpha=0.7,
               label="bag$\\to$stripping ($We_r{=}0.5\\sqrt{Re}$)")
    ax.axvline(check.WEBAG, ls=":", color="gray", lw=1)
    ax.set_xscale("log")
    ax.set_xlabel("$We_r = \\rho_g u^2 r/\\sigma$")
    ax.set_ylabel("initial $\\mathrm{d}d/\\mathrm{d}t$ [m/s]")
    ax.set_title(f"{CASE}: Reitz-Diwakar breakup rate vs [RD87] "
                 "(bag + stripping)")
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
