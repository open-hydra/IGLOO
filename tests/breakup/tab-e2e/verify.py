#!/usr/bin/env python3
"""Overlay IGLOO vs the ORA87 TAB first-breakup time (B-VAL-3, MOSE-style).

check.py is the GATE (analytic damped-oscillator onset + t_bu, xfail on A18); this
only DRAWS and writes OUTPUT/tab-e2e.svg. matplotlib imported defensively.

One panel: t_bu vs We_r (radius-based, ORA87 convention): the analytic crossing curve
(gated), the no-break region We_r<6 shaded, and IGLOO's measured brackets (error
bars). Post-A18-fix the points ride the analytic curve to sub-percent.
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
    ap = argparse.ArgumentParser(description="overlay IGLOO vs ORA87 t_bu(We_r)")
    ap.add_argument("--plot", action="store_true")
    args = ap.parse_args()
    if not os.path.isfile(TRAJ):
        print(f"[verify] {TRAJ} not found -- run the case first")
        return 0
    parts = check.load_trajectories(TRAJ)

    # analytic curve over the sweep range (slip from the design kv=0.5)
    slip = 0.5 * check.U_GAS
    wes, tbs = [], []
    for i in range(301):
        wer = 6.05 * (95.0 / 6.05) ** (i / 300.0)
        r = wer * check.SIGMA / (check.RHOG * slip**2)
        tb = check.tbu_analytic(r, slip)
        if tb is not None:
            wes.append(wer); tbs.append(tb * 1e3)

    meas = []
    for pid in sorted(parts):
        raw = sorted(parts[pid])          # record order is OMP-nondeterministic
        rows = [raw[0]]
        for r in raw[1:]:
            if r[0] > rows[-1][0]:
                rows.append(r)
        if len(rows) < 4:
            continue
        wer, tlo, thi, _ = check.drop_report(rows)
        if tlo is not None:
            meas.append((wer, 0.5 * (tlo + thi) * 1e3, 0.5 * (thi - tlo) * 1e3))

    fig, ax = plt.subplots(figsize=(6.2, 4.0))
    ax.plot(wes, tbs, "-", color="C0", lw=1.6,
            label="ORA87 damped oscillator: first $y{=}1$ crossing (gated)")
    if meas:
        w, t, e = zip(*meas)
        ax.errorbar(w, t, yerr=e, fmt="o", color="C1", ms=4, capsize=2,
                    label="IGLOO (measured bracket)")
    else:
        ax.text(0.45, 0.5, "no breakup events recorded",
                transform=ax.transAxes, ha="center", fontsize=9, color="C3")
    ax.axvspan(3, 6, color="gray", alpha=0.15, label="no-break: $We_r<6$ (gated)")
    ax.axvline(6, ls="--", color="gray", lw=0.8)
    ax.set_xscale("log")
    ax.set_xlim(3, 100)
    ax.set_xlabel("$We_r = \\rho_g u^2 r/\\sigma$"); ax.set_ylabel("$t_{bu}$ [ms]")
    ax.set_title(f"{CASE}: TAB first-breakup time vs ORA87 closed form [ORA87]")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(loc="upper left")
    fig.tight_layout()
    out = os.path.join(HERE, "OUTPUT", f"{CASE}.svg")
    fig.savefig(out, bbox_inches="tight", transparent=True)
    print(f"[verify] wrote {out}")
    if args.plot:
        plt.show()
    return 0


if __name__ == "__main__":
    sys.exit(main())
