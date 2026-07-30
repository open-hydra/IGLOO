#!/usr/bin/env python3
"""Overlay IGLOO vs the DIGITIZED MHB98 Fig. 2 M7 model line (water droplet, MOSE-style).

check.py is the GATE (IGLOO vs the LK kernel + IGLOO's beta vs the digitized M7 slope +
wet-bulb vs the M7 plateau); this only DRAWS OUTPUT/mhb98-water.svg. matplotlib is imported
defensively (no-op if absent).

IGLOO runs CEM+LK = the uniform-temperature Langmuir-Knudsen model = MHB98's M7 (M8 is the
finite-conductivity variant; the two are superimposed in Fig. 2, so only M7 was digitized).
Two panels vs the NORMALIZED D^2-law lifetime t/t_life, t_life = D0^2/K (K = LSQ slope,
extrapolated to D^2->0; NOT the last sample -- the water drop exits the box at d^2/d0^2~0.14):
  (a) d^2/d0^2 : IGLOO, its LK kernel (check.py's oracle), and the digitized M7 model line.
      Normalized (each by its own D0^2 and lifetime) because the paper's caption D0=1.1 mm
      (D0^2=1.21) disagrees with its own plot, whose M7 line starts at D0^2~1.10 (the axis even
      tops at 1.2). Since beta (Eq.28) is D0-INDEPENDENT and matches to ~1%, the two normalized
      d^2-law lines COINCIDE -- the honest presentation of an already-matching physics. beta is
      annotated in-panel (beta_IGLOO=6.44e-3 vs beta_M7=6.51e-3).
  (b) T_d [K]  : IGLOO wet-bulb history and the digitized M7 temperature line (same t/t_life).

The older eyeball Ranz-Marshall overlay (reference/mhb98_fig2.csv) is intentionally NOT
plotted: it reads systematically high (start 1.21, beta_exp~7.4e-3) whereas in the figure the
experimental points sit on M7 (~6.5e-3). A proper experiment re-digitization is a follow-up.
"""
import argparse
import os
import sys

import check                                    # reuse the gate's loaders + oracle

HERE = os.path.dirname(os.path.abspath(__file__))
CASE = os.path.basename(HERE)
TRAJ = os.path.join(HERE, "OUTPUT", "trajectories-A.dat")
REF_M7_D2 = os.path.join(HERE, "reference", "mhb98_fig2_M7_d2.csv")
REF_M7_TD = os.path.join(HERE, "reference", "mhb98_fig2_M7_Td.csv")

try:
    sys.path.insert(0, os.path.dirname(os.path.dirname(HERE)))
    import vv_style
    vv_style.apply()
    import matplotlib.pyplot as plt
except Exception as exc:
    print(f"[verify] plotting skipped ({exc.__class__.__name__}: {exc})")
    sys.exit(0)


def lsq_slope(ts, ys):
    """-slope of a straight-line fit y=a-K t (K>0 for a decaying d^2)."""
    n = len(ts); sx = sum(ts); sy = sum(ys)
    sxx = sum(t*t for t in ts); sxy = sum(t*y for t, y in zip(ts, ys))
    return -(n*sxy - sx*sy) / (n*sxx - sx*sx)


def igloo_and_kernel(parts):
    """(pid, t[s], d2_igloo[mm^2], d2_kernel[mm^2], Tp[K]) for the first sampled droplet."""
    for pid in sorted(parts):
        raw = sorted(parts[pid])          # record order is OMP-nondeterministic
        rows = [raw[0]]
        for r in raw[1:]:
            if r[0] > rows[-1][0]:
                rows.append(r)
        if len(rows) <= check.MIN_PTS:
            continue
        xs  = [r[0] for r in rows]
        Tps = [r[6] for r in rows]
        d2a = rows[0][7] ** 2                        # m^2
        pk = [d2a]
        for k in range(1, len(rows)):
            pk.append(check.rk4_span(xs[k-1], xs[k], Tps[k-1], Tps[k], pk[-1]))
        t   = [x / check.U_G for x in xs]
        MM2 = 1e6                                     # m^2 -> mm^2
        d2_ig = [(r[7] ** 2) * MM2 for r in rows]
        d2_kn = [max(v, 0.0) * MM2 for v in pk]
        return pid, t, d2_ig, d2_kn, Tps
    return None


def main():
    ap = argparse.ArgumentParser(description="overlay IGLOO vs MHB98 Fig.2 M7 (water)")
    ap.add_argument("--plot", action="store_true")
    args = ap.parse_args()
    if not os.path.isfile(TRAJ):
        print(f"[verify] {TRAJ} not found -- run the case first")
        return 0
    got = igloo_and_kernel(check.load_trajectories(TRAJ))
    if got is None:
        print("[verify] no sufficiently sampled droplet to plot")
        return 0
    pid, t, d2_ig, d2_kn, Tp = got
    tm7, ym7 = check.load_xy(REF_M7_D2)              # t[s], D^2[mm^2]
    ttd, ytd = check.load_xy(REF_M7_TD)              # t[s], Td[K]

    # beta (Eq.28, D0-independent) -- the quantitative match, annotated in panel (a).
    K_ig = lsq_slope(t, d2_ig)                       # mm^2/s (IGLOO absolute D^2 units)
    K_m7 = lsq_slope(tm7, ym7) if tm7 else float("nan")
    beta_ig = check.BETA_COEF * K_ig * 1e-6
    beta_m7 = check.BETA_COEF * K_m7 * 1e-6 if tm7 else float("nan")

    # Normalized D^2-law lifetime: t_life = D0^2 / K (extrapolate D^2->0), NOT the last
    # sample -- the water drop EXITS the box at d^2/d0^2~0.14, so t[-1] is not the lifetime.
    d0sq_ig = d2_ig[0]; life_ig = d0sq_ig / K_ig
    d0sq_m7 = ym7[0];   life_m7 = d0sq_m7 / K_m7 if tm7 else 1.0
    tn_ig = [tt / life_ig for tt in t]
    d2n_ig = [v / d0sq_ig for v in d2_ig]
    d2n_kn = [v / d0sq_ig for v in d2_kn]
    tn_m7 = [tt / life_m7 for tt in tm7] if tm7 else []
    d2n_m7 = [v / d0sq_m7 for v in ym7] if tm7 else []
    tn_td = [tt / life_m7 for tt in ttd] if ttd else []

    fig, ax = plt.subplots(1, 2, figsize=(9.8, 4.0))

    # Panel (a): normalized d^2/d0^2 vs t/t_life -- resolves the paper's own D0^2 caption
    # (1.1mm->1.21) vs plot (~1.10) inconsistency; equal beta => the two lines coincide.
    ax[0].plot(tn_ig, d2n_kn, "-", color="C0", lw=1.2, label="IGLOO LK kernel (gated oracle)")
    ax[0].plot(tn_ig, d2n_ig, "o", color="C1", ms=3.0, mfc="none", label="IGLOO")
    if tm7:
        ax[0].plot(tn_m7, d2n_m7, "s-", color="C2", ms=3.0, mfc="none", lw=1.0,
                   label="MHB98 M7 (digitized)")
    ax[0].set_xlabel(r"$t/t_\mathrm{life}$"); ax[0].set_ylabel("$d^2/d_0^2$")
    ax[0].set_title("water $d^2$-law (normalized)")
    ax[0].grid(True, alpha=0.25)
    rel_beta = abs(beta_ig - beta_m7) / beta_m7 * 100.0 if tm7 else float("nan")
    ax[0].text(0.04, 0.06,
               r"$\beta_{\mathrm{IGLOO}}=%.2f\times10^{-3}$" % (beta_ig*1e3) + "\n" +
               r"$\beta_{\mathrm{M7}}=%.2f\times10^{-3}$ (Eq.28, %.1f%%)"
               % (beta_m7*1e3, rel_beta),
               transform=ax[0].transAxes, va="bottom", ha="left")
    ax[0].legend(loc="upper right")

    ax[1].plot(tn_ig, Tp, "-", color="C3", lw=1.6, label="IGLOO")
    if ttd:
        ax[1].plot(tn_td, ytd, "s", color="C2", ms=3.5, mfc="none", label="MHB98 M7 (digitized)")
    ax[1].set_xlabel(r"$t/t_\mathrm{life}$"); ax[1].set_ylabel("$T_d$ [K]")
    ax[1].set_title("droplet temperature (wet-bulb)")
    ax[1].grid(True, alpha=0.25)
    ax[1].legend(loc="best")

    fig.suptitle(f"{CASE}: IGLOO vs MHB98 Fig.2 M7 [MHB98] (normalized; "
                 f"$T_G$=298 K, quiescent; ID {pid})")
    fig.tight_layout()
    out = os.path.join(HERE, "OUTPUT", f"{CASE}.svg")
    fig.savefig(out, bbox_inches="tight", transparent=True)
    print(f"[verify] wrote {out}  (beta_IGLOO={beta_ig:.3e}, beta_M7={beta_m7:.3e})")
    if args.plot:
        plt.show()
    return 0


if __name__ == "__main__":
    sys.exit(main())
