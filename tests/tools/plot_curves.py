#!/usr/bin/env python3
"""Render unit-family curve dumps into per-family overlay SVGs (MOSE-style).

The Fortran unit tests dump (x, production, reference...) curves via
support/verif_dump.f90 into <build_dir>/curves/<family>/<id>.dat with a
plot-metadata header (title/xlabel/ylabel/|-separated labels, optional
xlog/ylog). This script draws one figure per family (subplot grid, one panel
per dump) and writes unit-<family>.svg into --out-dir. Non-gating companion to
the e2e verify.py overlays; ctest runs the Fortran gates, not this file.
matplotlib is imported defensively (no-op if absent).

  ./plot_curves.py [--build-dir ../build/verif] [--out-dir ../build/verif/curves-svg]

Drawing convention (as the e2e overlays): data column 1 = production (open
orange markers), remaining columns = references (C0 solid, then C3/C2/C4
dashed).
"""
import argparse
import glob
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
TESTS = os.path.dirname(HERE)

try:
    sys.path.insert(0, TESTS)
    import vv_style                            # shared V&V style (CM/LaTeX math fonts)
    vv_style.apply()
    import matplotlib.pyplot as plt
except Exception as exc:
    print(f"[plot_curves] plotting skipped ({exc.__class__.__name__}: {exc})")
    sys.exit(0)

REF_STYLES = [("C0", "-"), ("C3", "--"), ("C2", "--"), ("C4", ":")]


def read_dump(path):
    """(meta dict, rows [[x, y1, ...], ...]) from one verif_dump file."""
    meta, rows = {}, []
    with open(path) as f:
        for line in f:
            s = line.strip()
            if not s:
                continue
            if s.startswith("#"):
                key, _, val = s.lstrip("# ").partition(":")
                meta[key.strip()] = val.strip()
                continue
            try:
                rows.append([float(v) for v in s.split()])
            except ValueError:
                continue
    return meta, rows


def draw_panel(ax, meta, rows):
    xs = [r[0] for r in rows]
    ncurves = min(len(r) for r in rows) - 1
    labels = [s.strip() for s in meta.get("labels", "").split("|")]
    labels += [f"curve {i}" for i in range(len(labels), ncurves)]
    for c in range(ncurves - 1, 0, -1):        # references first (under markers)
        color, ls = REF_STYLES[min(c - 1, len(REF_STYLES) - 1)]
        ax.plot(xs, [r[1 + c] for r in rows], ls, color=color, lw=1.4,
                label=labels[c], zorder=2)
    ax.plot(xs, [r[1] for r in rows], "o", color="C1", ms=3.0, mfc="none",
            label=labels[0], zorder=3)
    if meta.get("xlog") == "1":
        ax.set_xscale("log")
    if meta.get("ylog") == "1":
        ax.set_yscale("log")
    ax.set_xlabel(meta.get("xlabel", ""))
    ax.set_ylabel(meta.get("ylabel", ""))
    ax.set_title(meta.get("title", meta.get("id", "")))
    ax.grid(True, alpha=0.25)
    ax.legend(loc="best")


def main():
    ap = argparse.ArgumentParser(description="render unit-family curve dumps")
    ap.add_argument("--build-dir", default=None,
                    help="dir holding curves/ (default: build/verif/tests, the ctest cwd)")
    ap.add_argument("--out-dir", default=None)
    args = ap.parse_args()
    build_dir = args.build_dir or os.path.join(TESTS, "build", "verif", "tests")
    curves = os.path.join(build_dir, "curves")
    out_dir = args.out_dir or os.path.join(build_dir, "curves-svg")
    if not os.path.isdir(curves):
        print(f"[plot_curves] {curves} not found -- run ctest first")
        return 0
    os.makedirs(out_dir, exist_ok=True)

    n_fig = 0
    for fam in sorted(os.listdir(curves)):
        dumps = sorted(glob.glob(os.path.join(curves, fam, "*.dat")))
        if not dumps:
            continue
        n = len(dumps)
        ncol = min(3, n)
        nrow = math.ceil(n / ncol)
        fig, axs = plt.subplots(nrow, ncol, figsize=(4.6 * ncol, 3.4 * nrow),
                                squeeze=False)
        for k, path in enumerate(dumps):
            meta, rows = read_dump(path)
            if rows:
                draw_panel(axs[k // ncol][k % ncol], meta, rows)
        for k in range(n, nrow * ncol):        # blank the unused cells
            axs[k // ncol][k % ncol].set_axis_off()
        fig.suptitle(f"{fam} unit family: production vs literature", y=1.002)
        fig.tight_layout()
        out = os.path.join(out_dir, f"unit-{fam}.svg")
        fig.savefig(out, bbox_inches="tight", transparent=True)
        plt.close(fig)
        print(f"[plot_curves] wrote {out} ({n} panels)")
        n_fig += 1
    print(f"[plot_curves] {n_fig} family figures")
    return 0


if __name__ == "__main__":
    sys.exit(main())
