#!/usr/bin/env python3
"""tab-hf92 (etab-hf92 draws through this file): the products beside Hsiang & Faeth, drawn, never gated.

check.py gates the preconditions (P0-P3); this draws OUTPUT/<case>.svg from check.load() of the case directory:
  (a) [HF92] eq. [14] on its Fig. 10 axes, x = (rho_L/rho_G)^1/4 [mu_L/(rho_L d0 u0)]^1/2 We, y = rho_G SMD u0^2/sigma:
      the line y = 6.2 x with its +-20 % coefficient band, one marker per class (SMD = sum mdot / sum mdot/d);
  (b) SMD/d0 against We with eq. [14] and the regime onsets We = 13, 35, 80 ([HF93] p. 724);
  (c) MMD/SMD per class against the root-normal 1.2 ([HF92] p. 645, [HF93] Sec. 3.1);
  (d) TAB: the cumulative volume fraction against d/MMD of every class, pooled by regime, on [HF93] Fig. 1's axes
      (normal quantile against (d/MMD)^1/2) with the root-normal lines for MMD/SMD = 1.5, 1.2, 1.1 from the
      definition ((d/MMD)^1/2 normal by volume, cut at d/MMD = 0.01); ETAB: each class's sizes against its SMD;
  (e) [HF93] eq. [13], u0/u_b - 1 = 2.7 [(rho_G/rho_L)^1/2 d0/d]^2/3, per parcel at t_b = 5 t* (d and the slip u_b
      interpolated between the rows bracketing t_b), Fig. 6's data range shaded.
Parcels are taken by ascending ID, so the figure does not depend on the order of the records.
"""
import importlib.util
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
REGIMES = ((13.0, 35.0, "bag", "C0", "o"), (35.0, 80.0, "multimode", "C1", "s"), (80.0, math.inf, "shear", "C2", "^"))
RN_RATIOS = (1.5, 1.2, 1.1)          # [HF93] Fig. 1's legend order
RN_STYLES = ("--", "-", ":")
RN_XMIN = 0.1                        # (d/MMD)^1/2 cutoff: d/MMD >= 0.01
BAND = 0.2                           # [HF92] p. 649: 20 % on the coefficient of eq. [14]
C_UB, P_UB = 2.7, 2.0 / 3.0          # [HF93] eq. [13]
EQ13_FIT = (0.058, 1.07)             # [HF93] Fig. 6: the drawn correlation line
FIG6 = (0.067, 0.85)                 # [HF93] Fig. 6: the markers
P_TICKS = (0.1, 1, 5, 20, 50, 80, 95, 99)
D_TICKS = (0.0, 0.1, 0.5, 1.0, 2.0, 3.0, 4.0)
E_FRAME = ((1.0e-2, 1.0e1), (1.0e-2, 1.0e4))


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


check = load_module("tab_hf92_check", os.path.join(HERE, "check.py"))

try:
    sys.path.insert(0, os.path.dirname(os.path.dirname(HERE)))
    import vv_style
    vv_style.apply()
    import matplotlib.pyplot as plt
    from matplotlib.gridspec import GridSpec
except Exception as exc:
    print(f"[verify] plotting skipped ({exc.__class__.__name__}: {exc})")
    sys.exit(0)


def probit(p):
    """The standard normal quantile of p (bisection on erf)."""
    lo, hi = -9.0, 9.0
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if 0.5 * (1.0 + math.erf(mid / math.sqrt(2.0))) < p else (lo, mid)
    return 0.5 * (lo + hi)


def rn_ratio(s, n=2000):
    """MMD/SMD = E[1/x^2] for x = (d/MMD)^1/2 ~ N(1, s^2) by volume, x >= RN_XMIN (Simpson)."""
    a, b = RN_XMIN, 1.0 + 12.0 * s
    h = (b - a) / n
    num = den = 0.0
    for i in range(n + 1):
        x = a + i * h
        w = 1.0 if i in (0, n) else (4.0 if i % 2 else 2.0)
        g = math.exp(-0.5 * ((x - 1.0) / s) ** 2)
        num += w * g / (x * x)
        den += w * g
    return num / den


def rn_slope(ratio):
    lo, hi = 1.0e-3, 0.6
    for _ in range(40):
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if rn_ratio(mid) < ratio else (lo, mid)
    return 0.5 * (lo + hi)


def volume_positions(ds, ws):
    """(d, mid-step cumulative volume fraction) in ascending d, as check.mmd places them."""
    pts = sorted(zip(ds, ws))
    tot, acc, out = sum(w for _, w in pts), 0.0, []
    for d, w in pts:
        out.append((d, (acc + 0.5 * w) / tot))
        acc += w
    return out


def by_regime(ax, xs, ys, wes, label=True, ms=5.0):
    for lo, hi, name, col, mk in REGIMES:
        sel = [(x, y) for x, y, we in zip(xs, ys, wes) if lo <= we < hi]
        if sel:
            ax.plot(*zip(*sel), mk, color=col, ms=ms, mfc="none", mew=1.0,
                    label=f"IGLOO, {name}" if label else None)


def main(case_dir=HERE):
    name = os.path.basename(os.path.abspath(case_dir))
    if not os.path.isfile(os.path.join(case_dir, "OUTPUT", "outloc-A.dat")):
        print(f"[verify] {case_dir}/OUTPUT/outloc-A.dat not found -- run the case first")
        return 0
    case = check.load(case_dir)
    cl = [c for c in case["classes"] if c["ids"]]
    tab = case["model"].upper() == "TAB"
    rr, rg = case["rhol"] / case["rhog"], case["rhog"] / case["rhol"]
    wes = [c["we"] for c in cl]

    fig = plt.figure(figsize=(13.0, 8.6))
    gs = GridSpec(2, 6, figure=fig, hspace=0.36, wspace=1.1)
    ax_a, ax_b, ax_c = (fig.add_subplot(gs[0, 2 * i:2 * i + 2]) for i in range(3))
    ax_d, ax_e = fig.add_subplot(gs[1, 0:3]), fig.add_subplot(gs[1, 3:6])

    # (a) [HF92] eq. [14] on Fig. 10's axes
    xa = [rr ** 0.25 * math.sqrt(case["mul"] / (case["rhol"] * c["d0"] * c["u0"])) * c["we"] for c in cl]
    ya = [case["rhog"] * c["smd"] * c["u0"] ** 2 / case["sigma"] for c in cl]
    xl = [0.2, 20.0]
    ax_a.fill_between(xl, [(1 - BAND) * check.C_SMD * x for x in xl], [(1 + BAND) * check.C_SMD * x for x in xl],
                      color="0.85", lw=0, label=r"$\pm 20\,\%$ on the coefficient")
    ax_a.plot(xl, [check.C_SMD * x for x in xl], "k-", lw=1.2, label="[HF92] eq. [14]")
    by_regime(ax_a, xa, ya, wes)
    ax_a.set_xscale("log"); ax_a.set_yscale("log")
    ax_a.set_xlim(*xl); ax_a.set_ylim(0.5, 200.0)
    ax_a.set_xlabel(r"$(\rho_L/\rho_G)^{1/4}\,[\mu_L/(\rho_L d_0 u_0)]^{1/2}\,We$")
    ax_a.set_ylabel(r"$\rho_G\,\mathrm{SMD}\,u_0^2/\sigma$")
    ax_a.set_title("(a) SMD, [HF92] Fig. 10 axes")
    ax_a.grid(True, which="both", alpha=0.25)
    ax_a.legend(loc="upper left")

    # (b) SMD/d0 against We
    wl = [10.0 * (1000.0 / 10.0) ** (i / 200.0) for i in range(201)]
    d0, sg = cl[0]["d0"], case["sigma"]
    hf = [check.C_SMD * rr ** 0.25 * math.sqrt(case["mul"] / (case["rhol"] * d0 * math.sqrt(w * sg / (case["rhog"] * d0))))
          for w in wl]
    ax_b.fill_between(wl, [(1 - BAND) * v for v in hf], [(1 + BAND) * v for v in hf], color="0.85", lw=0)
    ax_b.plot(wl, hf, "k-", lw=1.2, label="[HF92] eq. [14]")
    by_regime(ax_b, wes, [c["smd"] / c["d0"] for c in cl], wes, label=False)
    for we, lab in check.ONSETS:
        ax_b.axvline(we, ls="--", color="gray", lw=0.8)
        ax_b.text(we * 1.06, 0.0062, lab, rotation=90, color="gray", fontsize=8, va="bottom")
    ax_b.set_xscale("log"); ax_b.set_yscale("log")
    ax_b.set_xlim(10.0, 1000.0); ax_b.set_ylim(0.005, 1.0)
    ax_b.set_xlabel(r"$We = \rho_G d_0 u_0^2/\sigma$")
    ax_b.set_ylabel(r"$\mathrm{SMD}/d_0$")
    ax_b.set_title(r"(b) SMD/$d_0$ against $We$")
    ax_b.grid(True, which="both", alpha=0.25)
    ax_b.legend(loc="upper right")

    # (c) MMD/SMD per class
    ax_c.axhline(1.2, color="k", lw=1.2, label="root-normal, MMD/SMD = 1.2")
    by_regime(ax_c, wes, [c["mmd"] / c["smd"] for c in cl], wes, label=False)
    ax_c.set_xscale("log")
    ax_c.set_xlim(10.0, 1000.0); ax_c.set_ylim(0.8, 1.6)
    ax_c.set_xlabel(r"$We$")
    ax_c.set_ylabel("MMD/SMD")
    ax_c.set_title("(c) MMD/SMD per class")
    ax_c.grid(True, which="both", alpha=0.25)
    ax_c.legend(loc="upper right")

    # (d) TAB: the size distribution on [HF93] Fig. 1's axes; ETAB: the class spread
    if tab:
        zl = (probit(0.0007), probit(0.995))
        for rat, ls in zip(RN_RATIOS, RN_STYLES):
            s = rn_slope(rat)
            ax_d.plot(zl, [1.0 + s * z for z in zl], ls, color="k", lw=1.1,
                      label=f"root-normal, MMD/SMD = {rat:.1f} ($s$ = {s:.3f})")
        for lo, hi, rname, col, mk in REGIMES:
            pts, ncl = [], 0
            for c in cl:
                if lo <= c["we"] < hi:
                    ncl += 1
                    ps = [case["parcels"][i] for i in c["ids"]]
                    pts += [(probit(v), math.sqrt(d / c["mmd"]))
                            for d, v in volume_positions([q["d"] for q in ps], [q["mdot"] for q in ps])]
            if pts:
                ax_d.plot(*zip(*pts), mk, color=col, ms=3.5, mfc="none", mew=0.8,
                          label=f"IGLOO, {rname} ({ncl} classes, {len(pts)} parcels)")
        ax_d.set_xlim(*zl); ax_d.set_ylim(0.0, 2.05)
        ax_d.set_xticks([probit(p / 100.0) for p in P_TICKS])
        ax_d.set_xticklabels([f"{p:g}" for p in P_TICKS])
        ax_d.set_yticks([math.sqrt(v) for v in D_TICKS])
        ax_d.set_yticklabels([f"{v:.1f}" for v in D_TICKS])
        ax_d.set_xlabel("cumulative volume fraction [%] (normal quantile scale)")
        ax_d.set_ylabel(r"$d/\mathrm{MMD}$ (square-root scale)")
        ax_d.set_title("(d) size distribution per class, [HF93] Fig. 1 axes")
        ax_d.legend(loc="upper left", fontsize=8)
    else:
        xs, ys, ws = [], [], []
        for c in cl:
            for i in c["ids"]:
                xs.append(c["we"]); ys.append(case["parcels"][i]["d"] / c["smd"]); ws.append(c["we"])
        by_regime(ax_d, xs, ys, ws, ms=4.0)
        spread = [(max(case["parcels"][i]["d"] for i in c["ids"]) - min(case["parcels"][i]["d"] for i in c["ids"]))
                  / c["smd"] for c in cl]
        k = max(range(len(cl)), key=lambda j: spread[j])
        ax_d.axhline(1.0, color="k", lw=0.8)
        ax_d.text(0.02, 0.04, f"(max - min)/SMD per class {min(spread):.2g} to {spread[k]:.3f} "
                  f"(We {cl[k]['we']:.0f})", transform=ax_d.transAxes, fontsize=9)
        ax_d.set_xscale("log")
        ax_d.set_xlim(10.0, 1000.0); ax_d.set_ylim(0.4, 1.6)
        ax_d.set_xlabel(r"$We$")
        ax_d.set_ylabel(r"$d/\mathrm{SMD}$ of the class")
        ax_d.set_title("(d) product sizes per class (no distribution law)")
        ax_d.grid(True, which="both", alpha=0.25)
        ax_d.legend(loc="upper left", fontsize=8)

    # (e) [HF93] eq. [13] at t_b = 5 t*
    xs, ys, ws, zero = [], [], [], 0
    for pid in sorted(case["parcels"]):
        q = case["parcels"][pid]
        c = case["classes"][q["cls"]]
        if q["tb"] is None or q["ub"] <= 0.0:
            zero += 1
            continue
        xs.append(math.sqrt(rg) * c["d0"] / q["tb"]["d"]); ys.append(c["u0"] / q["ub"] - 1.0); ws.append(c["we"])
    (x0, x1), (y0, y1) = E_FRAME
    out = zero + sum(1 for x, y in zip(xs, ys) if not (x0 <= x <= x1 and y0 <= y <= y1))
    ax_e.axvspan(*FIG6, color="0.88", lw=0, label="[HF93] Fig. 6 data range")
    xf = [EQ13_FIT[0] * (EQ13_FIT[1] / EQ13_FIT[0]) ** (i / 50.0) for i in range(51)]
    xo = [x0 * (x1 / x0) ** (i / 100.0) for i in range(101)]
    ax_e.plot(xo, [C_UB * x ** P_UB for x in xo], "k:", lw=0.9)
    ax_e.plot(xf, [C_UB * x ** P_UB for x in xf], "k-", lw=1.2, label="[HF93] eq. [13]")
    by_regime(ax_e, xs, ys, ws, ms=3.5)
    ax_e.set_xscale("log"); ax_e.set_yscale("log")
    ax_e.set_xlim(x0, x1); ax_e.set_ylim(y0, y1)
    ax_e.set_xlabel(r"$(\rho_G/\rho_L)^{1/2}\,d_0/d$ at $t_b = 5\,t^*$")
    ax_e.set_ylabel(r"$u_0/u_b - 1$")
    ax_e.set_title(r"(e) size and relative velocity at $t_b = 5\,t^*$")
    if out:
        ax_e.text(0.98, 0.04, f"{out} parcel(s) outside the frame ({zero} without a slip at $t_b$)",
                  transform=ax_e.transAxes, fontsize=8, ha="right")
    ax_e.grid(True, which="both", alpha=0.25)
    ax_e.legend(loc="upper left", fontsize=8)

    fig.suptitle(f"{name}: {case['model']} products of Hsiang & Faeth's water series "
                 f"($d_0$ = {cl[0]['d0'] * 1e3:g} mm, $We$ {min(wes):.0f}-{max(wes):.0f}) against [HF92] and [HF93]")
    out_svg = os.path.join(case_dir, "OUTPUT", f"{name}.svg")
    fig.savefig(out_svg, bbox_inches="tight", transparent=True)
    plt.close(fig)
    print(f"[verify] wrote {out_svg}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
