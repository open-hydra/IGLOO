#!/usr/bin/env python3
"""Plot IGLOO burn-box output vs the EXACT Beckstead reference (MOSE-style verify.py).

check.py is the GATE; this only draws the overlay and writes OUTPUT/burn-box.svg.
It never gates the suite. matplotlib is imported defensively (no-op if absent).

  ./verify.py [--plot]

Reference: reference/beckstead.txt (x, d[m], Tp[K]) -- the closed-form d^n(x) burn
(temperature-independent above T_ign) and the RK4 energy-balance Tp(x), both
reproduced from the declared law in check.py (not digitized). Two panels overlay
the IGLOO diameter (col 8) and temperature (col 7) as markers.
"""
import argparse
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
CASE = os.path.basename(HERE)
TRAJ = os.path.join(HERE, "OUTPUT", "trajectories-A.dat")
REF  = os.path.join(HERE, "reference", "beckstead.txt")
MIN_PTS = 5

try:
    import numpy as np
    sys.path.insert(0, os.path.dirname(os.path.dirname(HERE)))
    import vv_style                            # shared V&V style (CM/LaTeX math fonts)
    vv_style.apply()
    import matplotlib.pyplot as plt
except Exception as exc:
    print(f"[verify] plotting skipped ({exc.__class__.__name__}: {exc})")
    sys.exit(0)


def best_particle(path):
    """(xa, xs, d_um, Tp) of the best-sampled particle (monotone-x rows only)."""
    parts = {}
    with open(path) as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith("variables") or s.startswith("Zone"):
                continue
            c = s.split()
            if len(c) < 10:
                continue
            try:
                x, T, d = float(c[0]), float(c[6]), float(c[7])
                pid = int(c[9])
            except ValueError:
                continue
            parts.setdefault(pid, []).append((x, T, d))
    best = None
    for pid, raw in parts.items():
        raw.sort()
        rows = [raw[0]]
        for r in raw[1:]:
            if r[0] > rows[-1][0]:
                rows.append(r)
        if len(rows) >= MIN_PTS and (best is None or len(rows) > len(best)):
            best = rows
    if best is None:
        return None
    return (best[0][0], [r[0] for r in best],
            [r[2] * 1e6 for r in best], [r[1] for r in best])


def main():
    ap = argparse.ArgumentParser(description="overlay IGLOO vs Beckstead d^n + energy balance")
    ap.add_argument("--plot", action="store_true", help="also open a window")
    args = ap.parse_args()
    if not os.path.isfile(TRAJ):
        print(f"[verify] {TRAJ} not found -- run the case first")
        return 0
    ref = np.loadtxt(REF)
    xr, dr_um, Tr = ref[:, 0], ref[:, 1] * 1e6, ref[:, 2]
    got = best_particle(TRAJ)
    if got is None:
        print("[verify] no sufficiently sampled particle to plot")
        return 0
    xa, xs, d_um, Tp = got

    fig, axs = plt.subplots(2, 1, figsize=(7.0, 6.8))
    axs[0].plot(xr + xa, dr_um, "-", color="C0", lw=1.6, label="Beckstead:  $d^n(t)=d_0^n-K_{eff}\\,t$,  $K_{eff}=K_{burn}X_{eff}$", zorder=2)
    axs[0].plot(xs, d_um, "o", color="C1", ms=3.5, mfc="none", label="IGLOO", zorder=3)
    axs[0].set_ylabel("$d$ [$\\mu$m]")
    axs[0].set_title(f"{CASE}: Beckstead $d^n$ aluminium burn [Beck05]")
    axs[1].plot(xr + xa, Tr, "-", color="C0", lw=1.6, label="RK4:  $mc_p\\frac{dT_p}{dt}=\\pi d\\,\\mathrm{Nu}\\,k_g(T_g-T_p)+\\beta q_c|\\dot{m}|$,  $\\mathrm{Nu}=2$", zorder=2)
    axs[1].plot(xs, Tp, "o", color="C1", ms=3.5, mfc="none", label="IGLOO", zorder=3)
    axs[1].set_ylabel("$T_p$ [K]")
    axs[1].set_title("temperature (combustion heat release)")
    for ax in axs:
        ax.set_xlim(min(xs), max(xs))
        ax.set_xlabel("$x$ [m]")
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
