#!/usr/bin/env python3
"""Plot the swirl-wedge spirals: IGLOO parcels vs the exact cylindrical ODE.

check.py is the GATE; this only draws the overlay and writes OUTPUT/swirl-wedge.svg.
It never gates the suite. matplotlib is imported defensively (no-op if absent).

  ./verify.py [--plot]

Left: radius r(x) -- every parcel drifts outward (centrifugal, no pressure force on a
parcel); right: azimuthal velocity w_p(x) relaxing toward the local swirl omega*r.
Markers = IGLOO (after the fold: rows are in-sector); solid lines = the RK4 oracle
integrated from the exact injection state (reuses check.py's rhs and constants).
"""
import argparse
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
CASE = os.path.basename(HERE)
TRAJ = os.path.join(HERE, "OUTPUT", "trajectories-A.dat")

try:
    sys.path.insert(0, os.path.dirname(os.path.dirname(HERE)))
    import vv_style                            # shared V&V style (CM/LaTeX math)
    vv_style.apply()
    import matplotlib.pyplot as plt
    import check                               # reuse the gate oracle + constants
except Exception as exc:
    print(f"[verify] plotting skipped ({exc.__class__.__name__}: {exc})")
    sys.exit(0)

LABEL = {1: r"P1 co-rotating, $w_0=\Omega r_0$",
         2: r"P2 spin-up, $w_0=0$",
         3: r"P3 over-spun, $w_0=2\Omega r_0$"}


def main():
    ap = argparse.ArgumentParser(description="plot the swirl-wedge spirals")
    ap.add_argument("--plot", action="store_true", help="also open a window")
    args = ap.parse_args()
    if not os.path.isfile(TRAJ):
        print(f"[verify] {TRAJ} not found -- run the case first")
        return 0
    parts = check.load(TRAJ)

    fig, (ax_r, ax_w) = plt.subplots(1, 2, figsize=(9.2, 3.8))
    for i, pid in enumerate(sorted(parts)):
        rows = parts[pid]
        x0, y0, z0, u0, v0, w0 = rows[0][:6]
        r0 = math.hypot(y0, z0)
        s = (r0, 0.0, (y0 * v0 + z0 * w0) / r0, (y0 * w0 - z0 * v0) / r0)
        xs, rs, ws, ro, wo = [], [], [], [], []
        t_prev = 0.0
        for x, y, z, u, v, w, tp, d in rows:
            r = math.hypot(y, z)
            t = (x - check.X0) / check.U0
            s = check.rk4(check.cyl_rhs, s, t_prev, t, check.RK_H)
            t_prev = t
            xs.append(x); rs.append(r); ws.append((y * w - z * v) / r)
            ro.append(s[0]); wo.append(s[3])
        c = f"C{i % 10}"
        ax_r.plot(xs, ro, "-", color=c, lw=1.3, zorder=2)
        ax_r.plot(xs[::4], rs[::4], "o", color=c, ms=2.4, mfc="none", zorder=3,
                  label=LABEL.get(pid, f"ID {pid}"))
        ax_w.plot(xs, wo, "-", color=c, lw=1.3, zorder=2)
        ax_w.plot(xs[::4], ws[::4], "o", color=c, ms=2.4, mfc="none", zorder=3)
    ax_r.set_xlabel("$x$ [m]")
    ax_r.set_ylabel("$r$ [m]")
    ax_r.set_title(rf"{CASE}: solid-body swirl $\Omega={check.OMEGA}$, "
                   rf"$\mathrm{{St}}=\tau\Omega={check.OMEGA*check.TAU:.1f}$")
    ax_w.set_xlabel("$x$ [m]")
    ax_w.set_ylabel(r"$w_p$ [m/s]")
    ax_w.set_title(r"azimuthal velocity vs $\dot w=(\Omega r-w)/\tau-v_r w/r$")
    for ax in (ax_r, ax_w):
        ax.set_xlim(0, check.LX)
        ax.grid(True, alpha=0.25)
    ax_r.plot([], [], "-", color="0.3", label="cylindrical ODE (RK4)")
    ax_r.plot([], [], "o", color="0.3", mfc="none", label="IGLOO (2.5D wedge)")
    ax_r.legend(loc="center right", fontsize=8)
    fig.tight_layout()
    out = os.path.join(HERE, "OUTPUT", f"{CASE}.svg")
    fig.savefig(out, bbox_inches="tight", transparent=True)
    print(f"[verify] wrote {out}")
    if args.plot:
        plt.show()
    return 0


if __name__ == "__main__":
    sys.exit(main())
