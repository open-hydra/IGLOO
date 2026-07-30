#!/usr/bin/env python3
"""Overlay IGLOO vs PE87 T*(We) for the Pilch-Erdman breakup validation (MOSE-style).

check.py is the GATE (initial breakup rate vs the PE87 closed-form correlation at
We>=45, plus the d_stable plateau); this only DRAWS and writes
OUTPUT/pilch-erdman-e2e.svg -- it never gates the suite. matplotlib is imported
defensively (no-op if absent).

One panel: dimensionless total breakup time T* vs We (log-x, PE87 Fig. 7 layout):
  - the PE87 published correlation (5 regimes, continuous) -- the gated reference;
  - IGLOO points: each drop's measured initial rate inverted through the paper model
    (the T* that reproduces the measured dd/dt at the injection (d0, slip0)). All
    We>=45 points are gated post-A17-fix (2026-07-22).
No experimental shock-tube scatter is overlaid (PE87 Fig. 7 not digitized; the gate is
the closed form, so digitization would add provenance cost for a non-gating layer).
"""
import argparse
import os
import sys

import check                                    # reuse the gate's kernel + loader

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


def rate_at_tstar(d, slip, ts):
    """PE87 diameter-relaxation rate at (d, slip) with an IMPOSED T* (paper VdV/dStable)."""
    vdv = check.SQEPS * (0.75 * check.CD * ts + 3.0 * check.BB * ts**2)
    wec = 12.0 * (1.0 + 1.077 * (check.MUP / (check.RHOP * d * check.SIGMA) ** 0.5) ** 1.6)
    dstable = wec * check.SIGMA / (check.RHOG * slip**2 * (1.0 - vdv) ** 2 + 1e-30)
    tau = ts * d / (slip * check.SQEPS)
    return (dstable - d) / tau


def implied_tstar(d, slip, rate, lo=0.5, hi=60.0):
    """Invert rate -> T* by bisection (|rate| is monotone-decreasing in T*)."""
    if rate >= 0.0:
        return None
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        if rate_at_tstar(d, slip, mid) < rate:   # model still faster than measured
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def igloo_points(parts):
    """[(We0, T*_implied, gated?)] per drop, same selection rules as the gate."""
    pts = []
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
        we0 = d0 * check.RHOG * slip0**2 / check.SIGMA
        if we0 < 45.0:
            continue
        r_ig = check.measured_rate0(rows)
        if r_ig is None:
            continue
        ts = implied_tstar(d0, slip0, r_ig)
        if ts is not None:
            pts.append((we0, ts, we0 >= check.WE_MIN_GATE))
    return pts


def main():
    ap = argparse.ArgumentParser(description="overlay IGLOO vs PE87 T*(We)")
    ap.add_argument("--plot", action="store_true")
    args = ap.parse_args()
    if not os.path.isfile(TRAJ):
        print(f"[verify] {TRAJ} not found -- run the case first")
        return 0
    pts = igloo_points(check.load_trajectories(TRAJ))
    if not pts:
        print("[verify] no drops to plot")
        return 0

    wes = [13.0 * (3000.0 / 13.0) ** (i / 400.0) for i in range(401)]
    tss = [check.pe87_tstar(w) for w in wes]

    fig, ax = plt.subplots(figsize=(6.2, 4.0))
    ax.plot(wes, tss, "-", color="C0", lw=1.6,
            label="PE87 $T^*(We)$ correlation (gated)")
    g = [(w, t) for w, t, ok in pts if ok]
    i = [(w, t) for w, t, ok in pts if not ok]
    if g:
        ax.plot(*zip(*g), "o", color="C1", ms=4.5,
                label="IGLOO (implied $T^*$, gated)")
    if i:
        ax.plot(*zip(*i), "o", color="C1", ms=4.5, mfc="none",
                label="IGLOO (below gate scope)")
    ax.axvline(check.WE_MIN_GATE, ls=":", color="gray", lw=1)
    for wb in (18.0, 45.0, 351.0, 2670.0):
        ax.axvline(wb, ls="--", color="gray", lw=0.6, alpha=0.5)
    ax.set_xscale("log")
    ax.set_xlim(13, 3000); ax.set_ylim(0, 7)
    ax.set_xlabel("$We$"); ax.set_ylabel("$T^*$")
    ax.set_title(f"{CASE}: total breakup time vs PE87 Fig. 7 [PE87] "
                 "(25-drop sweep, slip$_0$=100 m/s)")
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
