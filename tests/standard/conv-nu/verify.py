#!/usr/bin/env python3
"""Plot IGLOO conv-nu output vs the EXACT Nu(Re) reference (MOSE-style verify.py).

check.py is the GATE; this only draws the overlay and writes OUTPUT/conv-nu.svg.
It never gates the suite. matplotlib is imported defensively (no-op if absent).

  ./verify.py [--plot]

Reference: reference/nu_re.txt (x, T_conv, T_cond). T_conv is the Ranz-Marshall
Nu(Re) closed form at the constant slip; T_cond is the Nu=2 conduction curve shown
dashed as the regression the case out-distinguishes. Both reproduced from the
declared law in check.py (not digitized). IGLOO T_p(x) (col 7) overlaid as markers.
"""
import argparse
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
CASE = os.path.basename(HERE)
TRAJ = os.path.join(HERE, "OUTPUT", "trajectories-A.dat")
REF  = os.path.join(HERE, "reference", "nu_re.txt")
VCOL, MIN_PTS = 6, 3

try:
    import numpy as np
    sys.path.insert(0, os.path.dirname(os.path.dirname(HERE)))
    import vv_style                            # shared V&V style (CM/LaTeX math fonts)
    vv_style.apply()
    import matplotlib.pyplot as plt
except Exception as exc:
    print(f"[verify] plotting skipped ({exc.__class__.__name__}: {exc})")
    sys.exit(0)


def best_particle(path, vcol):
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
                x, y = float(c[0]), float(c[vcol])
                pid = int(c[9])
            except ValueError:
                continue
            parts.setdefault(pid, []).append((x, y))
    best = max((sorted(r) for r in parts.values() if len(r) >= MIN_PTS),
               key=len, default=None)
    if best is None:
        return None
    return best[0][0], [r[0] for r in best], [r[1] for r in best]


def main():
    ap = argparse.ArgumentParser(description="overlay IGLOO vs Nu(Re) convective curve")
    ap.add_argument("--plot", action="store_true", help="also open a window")
    args = ap.parse_args()
    if not os.path.isfile(TRAJ):
        print(f"[verify] {TRAJ} not found -- run the case first")
        return 0
    ref = np.loadtxt(REF)
    xr, Tconv, Tcond = ref[:, 0], ref[:, 1], ref[:, 2]
    got = best_particle(TRAJ, VCOL)
    if got is None:
        print("[verify] no sufficiently sampled particle to plot")
        return 0
    xa, xs, Ts = got

    fig, ax = plt.subplots(figsize=(7.0, 4.2))
    ax.plot(xr + xa, Tconv, "-", color="C0", lw=1.6, label="Ranz-Marshall:  $\\mathrm{Nu}=2+0.6\\,\\mathrm{Re}^{1/2}\\,\\mathrm{Pr}^{1/3}$\n$T_p(x)=T_g+(T_a-T_g)\\,e^{-(x-x_a)/L}$,  $L=u_g\\tau_T$,  $\\tau_T=\\frac{\\rho_p c_p d^2}{6\\,\\mathrm{Nu}\\,k_g}$", zorder=2)
    ax.plot(xr + xa, Tcond, "--", color="C3", lw=1.2, label="$\\mathrm{Nu}=2$ conduction (regression)", zorder=1)
    ax.plot(xs, Ts, "o", color="C1", ms=3.5, mfc="none", label="IGLOO", zorder=3)
    ax.set_xlim(min(xs), max(xs))
    ax.set_xlabel("$x$ [m]")
    ax.set_ylabel("$T_p$ [K]")
    ax.set_title(f"{CASE}: Ranz-Marshall $\\mathrm{{Nu}}(\\mathrm{{Re}})$ at constant slip [RM52]")
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
