#!/usr/bin/env python3
"""Overlay IGLOO vs the digitized TC2012 Fig. 11 model curves (E-VAL-3, MOSE-style).

check.py is the GATE (variable-density TC mass rate + swelling); this only DRAWS and
writes OUTPUT/tc-hexadecane.svg. matplotlib imported defensively.

IGLOO's TC kernel IS TC2012's present model (eq. 16): its small-rate limit reduces to
eq. 16's stated form m_hat=(P_vs-chi)/[(T~s+1)/2] exactly (not the Stefan-Fuchs eq. 2b).
So this compares IGLOO to the RIGHT digitized curve.

Two panels vs a NORMALIZED lifetime t/t_life (each series scaled by its own final time):
  (a) d^2/d0^2 : IGLOO and the digitized TC2012 Fig. 11 present-model curve;
  (b) T/T_d0   : IGLOO wet-bulb history and the digitized TC2012 Fig. 11 curve.

Tier-1 property reconstruction (2026-07-28) -- see INFO.md/PROVENANCE.md:
  * Lv=2.58e5 = the effective latent heat of TC2012's OWN Table-1 psat curve (CC slope
    through its (T,Pvs/PN) anchors ~ Watson Lv at the ~490K wet-bulb); the earlier 2.9e5
    made CC-psat run 20-24% LOW vs Table 1, over-heating the drop (peak too early).
  * cp_l=2800 (properties.dat) = NIST/Chemeo n-hexadecane liquid Cp at the ~490K operating
    point; the earlier 2200 was the 298K value -- right property, wrong reference T.
With these, IGLOO reproduces the heating shape well (heat-frac ~0.54, matching Fig. 11) and
the plateau to ~0.7% (490.2 vs ~493.6 K).

The x-axis is NORMALIZED by each curve's lifetime because IGLOO's D_v (set by Le and the
case gas props) differs from TC2012's n-hexadecane D_v, so absolute tau=t*Dv/R^2 is not
comparable; the normalized SHAPE is the valid comparison. A residual remains in the d^2
mid-decline: a single scalar Lv cannot be BOTH the psat-curve slope (~258 kJ/kg) and the
optimal energy sink (~227), so IGLOO evaporates a touch fast. Closing it needs a psat curve
DECOUPLED from the sink Lv (Tier 2, deferred -- see plan-bucket/tc-hexadecane-tier2-*.md).
Both digitized curves are NON-gating; the tight gate (check.py) is IGLOO-vs-the-TC-kernel.
"""
import argparse
import os
import sys

import check

HERE = os.path.dirname(os.path.abspath(__file__))
CASE = os.path.basename(HERE)
TRAJ = os.path.join(HERE, "OUTPUT", "trajectories-A.dat")
REF_D2 = os.path.join(HERE, "reference", "tc2012_fig11_d2.csv")
REF_T  = os.path.join(HERE, "reference", "tc2012_fig11_T.csv")

R_D0  = 0.5 * check.D0                       # initial drop radius = 10 um
TAU_PER_X = check.DV / (check.U_G * R_D0**2) # tau = t*Dv/R_d0^2, t = x/u_g

try:
    sys.path.insert(0, os.path.dirname(os.path.dirname(HERE)))
    import vv_style
    vv_style.apply()
    import matplotlib.pyplot as plt
except Exception as exc:
    print(f"[verify] plotting skipped ({exc.__class__.__name__}: {exc})")
    sys.exit(0)


def load_csv(path):
    xs, ys = [], []
    if os.path.isfile(path):
        for line in open(path):
            line = line.strip()
            if not line:
                continue
            a, b = line.split()
            xs.append(float(a)); ys.append(float(b))
    return xs, ys


def main():
    ap = argparse.ArgumentParser(description="overlay IGLOO vs TC2012 Fig.11")
    ap.add_argument("--plot", action="store_true")
    args = ap.parse_args()
    if not os.path.isfile(TRAJ):
        print(f"[verify] {TRAJ} not found -- run the case first")
        return 0
    parts = check.load_trajectories(TRAJ)
    pid = next((p for p in sorted(parts) if len(parts[p]) >= 6), None)
    if pid is None:
        print("[verify] no sufficiently sampled droplet")
        return 0
    rows = sorted(parts[pid])
    seen = [rows[0]]
    for r in rows[1:]:
        if r[0] > seen[-1][0]:
            seen.append(r)
    d0 = seen[0][3]
    tmax = seen[-1][0] or 1.0
    tn  = [r[0] / tmax for r in seen]            # IGLOO normalized lifetime
    d2n = [(r[3] / d0)**2 for r in seen]
    Tn  = [r[2] / check.TP0 for r in seen]

    tau_d2, ref_d2 = load_csv(REF_D2)
    tau_T,  ref_T  = load_csv(REF_T)
    n_d2 = (max(tau_d2) if tau_d2 else 1.0) or 1.0
    n_T  = (max(tau_T) if tau_T else 1.0) or 1.0

    fig, ax = plt.subplots(1, 2, figsize=(9.4, 3.8))
    ax[0].plot(tn, d2n, "-", color="C1", lw=1.6, label="IGLOO")
    if tau_d2:
        ax[0].plot([t / n_d2 for t in tau_d2], ref_d2, "s", color="k", ms=3.5, mfc="none",
                   label="TC2012 Fig. 11 (digitized)")
    ax[0].axhline(1.0, ls=":", color="gray", lw=0.8)
    ax[0].set_xlabel(r"$t/t_\mathrm{life}$ (normalized)"); ax[0].set_ylabel("$d^2/d_0^2$")
    ax[0].set_title("drop size (swelling then $D^2$-law)"); ax[0].grid(True, alpha=0.25)
    ax[0].legend(loc="best")

    ax[1].plot(tn, Tn, "-", color="C3", lw=1.6, label="IGLOO")
    if tau_T:
        ax[1].plot([t / n_T for t in tau_T], ref_T, "s", color="k", ms=3.5, mfc="none",
                   label="TC2012 Fig. 11 (digitized)")
    ax[1].set_xlabel(r"$t/t_\mathrm{life}$ (normalized)"); ax[1].set_ylabel("$T/T_{d,0}$")
    ax[1].set_title("drop temperature (wet-bulb plateau)"); ax[1].grid(True, alpha=0.25)
    ax[1].legend(loc="best")

    fig.suptitle(f"{CASE}: n-hexadecane TC evaporation vs TC2012 Fig. 11 "
                 "($R_{d,0}$=10 um, $T_\\infty$=600 K)")
    fig.tight_layout()
    out = os.path.join(HERE, "OUTPUT", f"{CASE}.svg")
    fig.savefig(out, bbox_inches="tight", transparent=True)
    print(f"[verify] wrote {out}")
    if args.plot:
        plt.show()
    return 0


if __name__ == "__main__":
    sys.exit(main())
