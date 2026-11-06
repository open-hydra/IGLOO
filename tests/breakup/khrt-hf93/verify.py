#!/usr/bin/env python3
"""khrt-hf93: Reitz KH-RT's parent drop beside Hsiang & Faeth's core drop, drawn, never gated (OUTPUT/khrt-hf93.svg).

check.py gates the preconditions (H0, H1); this draws from check.load(), one parcel per class (its lowest ID: a class's
parcels are copies of one deterministic path):
  (a) d/d0 against We: the parent before its first RT (the d of the last record, a trajectory row or a plane record,
      before the RT row), the parcel at t_b = 5 t* (hollow, interpolated between the rows that bracket t_b), and
      [HF93] d_b/d0 = (1 + C)(4 Eo_cr/(3 C_D We))^1/2 with the coefficient 2.0656 its constants give and with eq. [8]'s
      printed 2.3;
  (b) the relative velocity over u0 at the first RT and at t_b against [HF93] u_b/u0 = 1/(1 + C);
  (c) t/t* of the first RT against [HF93]'s t_b/t* = 5.
C = (3 C_D t_b/(4 t*))(rho_G/rho_L)^1/2 ([HF93] eq. [4]), C_D = 5, Eo_cr = 16, t_b/t* = 5 (p. 726), the densities of
the case. The correlations are drawn solid over the shear-breakup range We 100-600 their fit covers and dotted below
it; the regime onsets We = 13, 35, 80 ([HF93] p. 724) are marked.
"""
import importlib.util
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
CD, EO_CR, TB_STAR = 5.0, 16.0, 5.0          # [HF93] p. 726
K_PRINTED = 2.3                              # [HF93] eq. [8]'s printed coefficient
WE_FIT = (100.0, 600.0)                      # [HF93] p. 726: the shear breakup regime of the fit
REGIMES = ((13.0, 35.0, "bag", "C0", "o"), (35.0, 80.0, "multimode", "C1", "s"), (80.0, math.inf, "shear", "C2", "^"))


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


check = load_module("khrt_hf93_check", os.path.join(HERE, "check.py"))

try:
    sys.path.insert(0, os.path.dirname(os.path.dirname(HERE)))
    import vv_style
    vv_style.apply()
    import matplotlib.pyplot as plt
except Exception as exc:
    print(f"[verify] plotting skipped ({exc.__class__.__name__}: {exc})")
    sys.exit(0)


def by_regime(ax, xs, ys, wes, filled, label):
    for lo, hi, name, col, mk in REGIMES:
        sel = [(x, y) for x, y, we in zip(xs, ys, wes) if lo <= we < hi]
        if sel:
            ax.plot(*zip(*sel), mk, color=col, ms=5.0, mfc=col if filled else "none", mew=1.0,
                    label=f"{label}, {name}" if label else None)


def ranged(ax, f, style, color, label, we_lo):
    """f(We) solid (or dashed) over the fit's range, dotted from we_lo up to it."""
    lo = [we_lo * (WE_FIT[0] / we_lo) ** (i / 60.0) for i in range(61)]
    hi = [WE_FIT[0] * (WE_FIT[1] / WE_FIT[0]) ** (i / 60.0) for i in range(61)]
    ax.plot(lo, [f(w) for w in lo], ":", color=color, lw=1.0)
    ax.plot(hi, [f(w) for w in hi], style, color=color, lw=1.3, label=label)


def onsets(ax, y):
    for we, lab in check.hf92.ONSETS:
        ax.axvline(we, ls="--", color="gray", lw=0.8)
        ax.text(we * 1.06, y, lab, rotation=90, color="gray", fontsize=8, va="bottom")


def main(case_dir=HERE):
    name = os.path.basename(os.path.abspath(case_dir))
    if not os.path.isfile(os.path.join(case_dir, "OUTPUT", "planes-A.dat")):
        print(f"[verify] {case_dir}/OUTPUT/planes-A.dat not found -- run the case first")
        return 0
    case = check.load(case_dir)
    rg = case["rhog"] / case["rhol"]
    c_hf = 0.75 * CD * TB_STAR * math.sqrt(rg)
    k_hf = math.sqrt(4.0 * EO_CR / (3.0 * CD))
    cl = [c for c in case["classes"] if c["ids"] and case["parcels"][c["ids"][0]]["rt"]]
    qs = [case["parcels"][c["ids"][0]] for c in cl]
    wes = [c["we"] for c in cl]
    tb = [(c["we"], q) for c, q in zip(cl, qs) if q["tb"] is not None]
    we_lo = 10.0

    fig, (ax_a, ax_b, ax_c) = plt.subplots(1, 3, figsize=(15.0, 4.8))
    fig.subplots_adjust(wspace=0.3)

    # (a) the drop size
    ranged(ax_a, lambda w: (1.0 + c_hf) * k_hf / math.sqrt(w), "-", "k",
           rf"[HF93] eqs. [7]-[8], $(4Eo_{{cr}}/3\bar C_D)^{{1/2}}$ = {k_hf:.4f}", we_lo)
    ranged(ax_a, lambda w: (1.0 + c_hf) * K_PRINTED / math.sqrt(w), "--", "0.45",
           rf"[HF93] eq. [8] as printed, {K_PRINTED}", we_lo)
    by_regime(ax_a, wes, [q["rt"]["d_before"] / c["d0"] for c, q in zip(cl, qs)], wes, True, "IGLOO before its first RT")
    by_regime(ax_a, [w for w, _ in tb], [q["tb"]["d"] / c["d0"] for c, (_, q) in zip(cl, tb)],
              [w for w, _ in tb], False, "IGLOO at $t_b = 5\\,t^*$")
    onsets(ax_a, 1.2e-3)
    ax_a.set_xscale("log"); ax_a.set_yscale("log")
    ax_a.set_xlim(we_lo, 1000.0); ax_a.set_ylim(1.0e-3, 2.0)
    ax_a.set_xlabel(r"$We = \rho_G d_0 u_0^2/\sigma$")
    ax_a.set_ylabel(r"$d/d_0$")
    ax_a.set_title(r"(a) the parent drop against the core drop $d_b$")
    ax_a.grid(True, which="both", alpha=0.25)
    ax_a.legend(loc="upper center", bbox_to_anchor=(0.5, -0.16), fontsize=8)

    # (b) the relative velocity
    ranged(ax_b, lambda w: 1.0 / (1.0 + c_hf), "-", "k", rf"[HF93] eqs. [2]-[4], $1/(1 + C)$ = {1.0 / (1.0 + c_hf):.4f}",
           we_lo)
    by_regime(ax_b, wes, [q["rt"]["slip"] / c["u0"] for c, q in zip(cl, qs)], wes, True, "IGLOO at its first RT")
    by_regime(ax_b, [w for w, _ in tb], [q["ub"] / c["u0"] for c, (_, q) in zip(cl, tb)], [w for w, _ in tb], False,
              "IGLOO at $t_b = 5\\,t^*$")
    onsets(ax_b, 0.04)
    ax_b.set_xscale("log")
    ax_b.set_xlim(we_lo, 1000.0); ax_b.set_ylim(0.0, 1.1)
    ax_b.set_xlabel(r"$We$")
    ax_b.set_ylabel(r"$|u_g - u_p|/u_0$")
    ax_b.set_title("(b) the relative velocity")
    ax_b.grid(True, which="both", alpha=0.25)
    ax_b.legend(loc="upper center", bbox_to_anchor=(0.5, -0.16), fontsize=8)

    # (c) the time of the first RT
    ranged(ax_c, lambda w: TB_STAR, "-", "k", r"[HF93] $t_b/t^* = 5$", we_lo)
    by_regime(ax_c, wes, [q["rt"]["t"] / c["tstar"] for c, q in zip(cl, qs)], wes, True, "IGLOO, first RT")
    onsets(ax_c, 0.15)
    ax_c.set_xscale("log")
    ax_c.set_xlim(we_lo, 1000.0); ax_c.set_ylim(0.0, 6.0)
    ax_c.set_xlabel(r"$We$")
    ax_c.set_ylabel(r"$t/t^*$, $t^* = d_0 (\rho_L/\rho_G)^{1/2}/u_0$")
    ax_c.set_title("(c) the time of the first RT")
    ax_c.grid(True, which="both", alpha=0.25)
    ax_c.legend(loc="upper center", bbox_to_anchor=(0.5, -0.16), fontsize=8)

    fig.suptitle(f"{name}: Reitz KH-RT's parent drop on Hsiang & Faeth's water series ($d_0$ = {cl[0]['d0'] * 1e3:g} mm, "
                 f"$We$ {min(wes):.0f}-{max(wes):.0f}) against [HF93]'s core drop at the end of breakup")
    out_svg = os.path.join(case_dir, "OUTPUT", f"{name}.svg")
    fig.savefig(out_svg, bbox_inches="tight", transparent=True)
    plt.close(fig)
    print(f"[verify] wrote {out_svg}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
