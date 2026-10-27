#!/usr/bin/env python3
"""Plot IGLOO's burnout times beside Beckstead's Fig. 1 data and Eq. 7 (MOSE-style verify.py).

check.py is the GATE; this only draws OUTPUT/beckstead-air.svg and never gates the suite. matplotlib is imported
defensively (no-op if absent).

  ./verify.py [--plot]

On Fig. 1's axes (log-log, t_b [ms] against D [um], 10-1000 um x 0.1-10^4 ms): the 183 markers digitized from
[Beck05] Fig. 1 (reference/beckstead_fig1.csv, one grey glyph: the source labels of the stacked markers are not
recovered), the figure's fit as drawn (t = 0.00258 D^1.989, dashed; its label prints the rounded 0.003 D^1.99),
Eq. 7 at air, 1 atm, 300 K (check.py's K_eff, solid), and IGLOO's 25 burnout times (the exit records' t, markers) at
their injection diameters (bc.txt 2 rp). The parcels are taken by ascending ID, never in file order.
"""
import argparse
import os
import sys

import check                                    # the gate's Eq. 7 constants and readers

HERE = os.path.dirname(os.path.abspath(__file__))
CASE = os.path.basename(HERE)
REF = os.path.join(HERE, "reference", "beckstead_fig1.csv")
FIT_A, FIT_N = 0.00258, 1.989                     # Fig. 1's fit as drawn (the label prints 0.003 D^1.99)

try:
    sys.path.insert(0, os.path.dirname(os.path.dirname(HERE)))
    import vv_style                            # shared V&V style (CM/LaTeX math fonts)
    vv_style.apply()
    import matplotlib.pyplot as plt
except Exception as exc:
    print(f"[verify] plotting skipped ({exc.__class__.__name__}: {exc})")
    sys.exit(0)


def fig1_points(path):
    """(D [um], t [ms]) of every digitized marker; '#' lines skipped."""
    ds, ts = [], []
    for line in open(path):
        s = line.strip()
        if not s or s.startswith("#"):
            continue
        c = [x.strip() for x in s.split(",")]
        ds.append(float(c[0])); ts.append(float(c[1]))
    return ds, ts


def igloo_points():
    """(ID, D0 [um], t_exit [ms]) of every exit record, by ascending ID."""
    d0 = check.inlet_d0(check.BC)
    out = []
    for e in check.rows(check.OUTLOC):
        c = check.cell_of(e["Y"], e["Z"])
        if c is not None:
            out.append((e["ID"], d0[c] * 1e6, e["t"] * 1e3))
    return sorted(out)


def main():
    ap = argparse.ArgumentParser(description="overlay IGLOO burnout times vs Beckstead Fig. 1 and Eq. 7")
    ap.add_argument("--plot", action="store_true", help="also open a window")
    args = ap.parse_args()
    if not os.path.isfile(check.OUTLOC):
        print(f"[verify] {check.OUTLOC} not found -- run the case first")
        return 0
    pts = igloo_points()
    if not pts:
        print("[verify] no exit record to plot")
        return 0
    dexp, texp = fig1_points(REF)
    keff_ms = check.KEFF / check.UNIT                # um^1.8/ms
    dl = [10.0 * 100.0 ** (k / 200.0) for k in range(201)]

    fig, ax = plt.subplots(figsize=(6.8, 5.6))
    ax.loglog(dexp, texp, "o", color="0.6", ms=4.0, mfc="0.75", mew=0.6,
              label="Beckstead (2005) Fig. 1, 12 sources (digitized)", zorder=1)
    ax.loglog(dl, [FIT_A * d ** FIT_N for d in dl], "--", color="k", lw=1.2,
              label=r"Fig. 1 fit as drawn, $t=0.00258\,D^{1.989}$ (label: $0.003\,D^{1.99}$)", zorder=2)
    ax.loglog(dl, [d ** check.NB / keff_ms for d in dl], "-", color="C0", lw=1.6,
              label=r"Eq. 7, air, 1 atm, 300 K: $t_b=a\,D^{1.8}/(X_\mathrm{eff}\,p^{0.1}\,T_0^{0.2})$", zorder=3)
    ax.loglog([p[1] for p in pts], [p[2] for p in pts], "s", color="C3", ms=5.0, mfc="none", mew=1.2,
              label="IGLOO burnout times (exit records)", zorder=4)
    ax.set_xlim(10.0, 1000.0)
    ax.set_ylim(0.1, 1.0e4)
    ax.set_xlabel(r"$D$ [$\mu$m]")
    ax.set_ylabel(r"$t_b$ [ms]")
    ax.set_title(f"{CASE}: aluminium burning times [Beck05]")
    ax.grid(True, which="major", alpha=0.25)
    ax.legend(loc="upper left")
    fig.tight_layout()
    out = os.path.join(HERE, "OUTPUT", f"{CASE}.svg")
    fig.savefig(out, bbox_inches="tight", transparent=True)
    lo, hi = pts[0], pts[-1]
    print(f"[verify] wrote {out}  ({len(pts)} IGLOO points, {len(dexp)} digitized; "
          f"ID {lo[0]}: D0 {lo[1]:.3f} um, {lo[2]:.4f} ms; ID {hi[0]}: D0 {hi[1]:.3f} um, {hi[2]:.1f} ms)")
    if args.plot:
        plt.show()
    return 0


if __name__ == "__main__":
    sys.exit(main())
