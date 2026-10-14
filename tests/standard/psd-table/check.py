#!/usr/bin/env python3
"""psd-table / psd-normal / psd-lognormal / psd-rosinrammler: the injected diameters follow the
size distribution of bc.txt column 8 (this file is shared through symlinks; the case is the law
named in its own INPUT/bc.txt).

The box (1 x 200 x 100 cells) injects in the cell-centre mode (fsample = 1, no ds, no mdot-max): one
parcel per inlet cell, 20 000 parcels, each carrying its cell's mass flow -- equal on the uniform
inlet, so the unweighted distribution of the parcel diameters is the injected mass-flow distribution.
Each parcel's diameter and position come from its FIRST trajectory row (written at injection), its
mass flow from its exit record, joined by ID. The prescribed CDFs are typed below from the case
inputs, never read back from the solver.

Gates (alpha = 1e-3 for both statistics, ks_ad.py):
  KS  Kolmogorov-Smirnov D_n < 1.94947/sqrt(n)        AD  Anderson-Darling A^2 < 5.969 (case 0)
  N   20 000 parcels, IDs 1..20000 once each in both files, one first row per inlet cell centre at x = 0
  M   every parcel's mass flow equal (relative spread <= 1e-12); the mass-weighted D_n equal to D_n
  L   bc.txt columns 6-8 (and table6.dat) equal to the typed constants: the oracle and the fixture agree
  W   the cell-centre mode ran (witness line) and, for the table, one notice that column 6 lies outside
      the table (the table, not rp, sets the size)
"""
import math
import os
import sys
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "tools"))
from ks_ad import K_CRIT, A2_CRIT, ks_d, ks_d_weighted, ad_a2, ks_pvalue, ad_pvalue  # noqa: E402

TRAJ, EXIT, BC, LOG = "OUTPUT/trajectories-A.dat", "OUTPUT/outloc-A.dat", "INPUT/bc.txt", "run_out.txt"
ALPHA = 1.0e-3
NY, NZ, LY, LZ = 200, 100, 0.05, 0.05
NPARCEL = NY * NZ

# ---- the laws, typed from the case inputs (bc.txt col 6 = rp, col 7 = sigmap, col 8 = law) ----
R_TAB = [2.5e-6, 4.0e-6, 6.0e-6, 1.5e-5, 2.25e-5, 3.0e-5]     # table6.dat column 1, radius [m]
M_TAB = [0.05, 0.10, 0.15, 0.35, 0.25, 0.10]                   # table6.dat column 2, mass fraction


def phi(z):
    return 0.5 * math.erfc(-z / math.sqrt(2.0))


def table_cdf():
    """Histogram reading: class k uniform between the mid-points to its neighbours, the end classes
    extended by half the adjacent spacing; edges on the radius, the CDF on the diameter 2r."""
    K = len(R_TAB)
    e = [R_TAB[0] - 0.5 * (R_TAB[1] - R_TAB[0])]
    e += [0.5 * (R_TAB[k] + R_TAB[k + 1]) for k in range(K - 1)]
    e += [R_TAB[-1] + 0.5 * (R_TAB[-1] - R_TAB[-2])]
    e = [2.0 * x for x in e]
    tot = sum(M_TAB)
    c = [0.0]
    for m in M_TAB:
        c.append(c[-1] + m / tot)

    def F(d):
        if d <= e[0]:
            return 0.0
        if d >= e[-1]:
            return 1.0
        k = max(i for i in range(K) if e[i] <= d)
        return c[k] + (d - e[k]) / (e[k + 1] - e[k]) * (c[k + 1] - c[k])
    return F


def normal_cdf(mu, sig):
    """Normal truncated to d > 0 (the sampler rejects d <= 0)."""
    p0 = phi(-mu / sig)
    return lambda d: 0.0 if d <= 0.0 else (phi((d - mu) / sig) - p0) / (1.0 - p0)


def lognormal_cdf(mean, std):
    s2 = math.log(1.0 + (std / mean) ** 2)
    m, s = math.log(mean) - 0.5 * s2, math.sqrt(s2)
    return lambda d: 0.0 if d <= 0.0 else phi((math.log(d) - m) / s)


def rosinrammler_cdf(x0, n):
    return lambda d: 0.0 if d <= 0.0 else 1.0 - math.exp(-(d / x0) ** n)


# token -> (rp, sigmap, description, CDF)
LAWS = {
    "INPUT/table6.dat": (1.0e-4, 0.0, "table6.dat (radius 2.5-30 um, bimodal)", table_cdf()),
    "Normal":           (1.0e-5, 1.0e-5, "Normal mu = 20 um, sigma = 10 um, truncated to d > 0",
                         normal_cdf(2.0e-5, 1.0e-5)),
    "LogNormal":        (1.0e-5, 1.0e-5, "LogNormal mean 20 um, std 10 um", lognormal_cdf(2.0e-5, 1.0e-5)),
    "RosinRammler":     (1.5e-5, 2.5, "Rosin-Rammler x0 = 30 um, n = 2.5", rosinrammler_cdf(3.0e-5, 2.5)),
}


def rows(path):
    for line in open(path):
        try:
            yield [float(v) for v in line.split()]
        except ValueError:
            continue


def read_table6(path):
    out = []
    for line in open(path):
        s = line.strip()
        if s and s[0] not in "#!;":
            out.append(tuple(float(v) for v in s.split()))
    return out


def main():
    for f in (TRAJ, EXIT, BC, LOG):
        if not os.path.exists(f):
            print(f"[FAIL] {f} missing")
            return 1
    log = open(LOG, errors="replace").read()
    payload = open(BC).readlines()[1].split()                  # first record: face 1, an inlet
    token = payload[7]
    if token not in LAWS:
        print(f"[FAIL] bc.txt column 8 = {token!r}: no oracle typed for it")
        return 1
    rp, sigmap, desc, F = LAWS[token]

    first = {}
    for r in rows(TRAJ):
        if len(r) >= 10:
            pid = int(r[9])
            if pid not in first or r[0] < first[pid][0]:
                first[pid] = r
    mdot = {}
    for r in rows(EXIT):
        if len(r) >= 9:
            mdot[int(r[8])] = r[6]
    ids = sorted(first)
    d = [first[i][7] for i in ids]
    cells = {(round(first[i][1] / (LY / NY) - 0.5), round(first[i][2] / (LZ / NZ) - 0.5)) for i in ids}
    w = [mdot.get(i, float("nan")) for i in ids]

    n = len(d)
    Dn, A2 = ks_d(d, F), ad_a2(d, F)
    Dcrit, A2crit = K_CRIT[ALPHA] / math.sqrt(n), A2_CRIT[ALPHA]
    wf = [Fraction(x) for x in w] if all(math.isfinite(x) for x in w) else None
    Dw = ks_d_weighted(d, [float(x) for x in wf], F) if wf else float("nan")
    spread = (max(w) - min(w)) / (sum(w) / n) if wf else float("inf")

    checks = [
        ("W  cell-centre injection ran (fsample = 1, one parcel per inlet cell)",
         "Placing particles at the bounday every  1 cells" in log and f"number of particles = {NPARCEL:6d}" in log),
        (f"N  {n} parcels, IDs 1..{NPARCEL} once each in the trajectory and exit files",
         n == NPARCEL and ids == list(range(1, NPARCEL + 1)) and sorted(mdot) == ids),
        (f"N  first rows at x = 0, one per inlet cell centre ({len(cells)} distinct of {NPARCEL})",
         len(cells) == NPARCEL and all(first[i][0] == 0.0 for i in ids)
         and all(0 <= j < NY and 0 <= k < NZ for j, k in cells)),
        (f"M  every parcel's mass flow equal: relative spread {spread:.3e} <= 1e-12", spread <= 1.0e-12),
        (f"M  mass-weighted D_n {Dw:.6e} equal to D_n {Dn:.6e}", abs(Dw - Dn) <= 1.0e-12),
        (f"L  bc.txt columns 6-8 = {rp:.5E} {sigmap:.5E} {token}",
         float(payload[5]) == rp and float(payload[6]) == sigmap),
        (f"KS {desc}: D_n = {Dn:.4e} < {Dcrit:.4e} (alpha {ALPHA:g}; p = {ks_pvalue(Dn, n):.3f})", Dn < Dcrit),
        (f"AD {desc}: A^2 = {A2:.4f} < {A2crit} (alpha {ALPHA:g}; p = {ad_pvalue(A2):.3f})", A2 < A2crit),
    ]
    if token.endswith("table6.dat"):
        tab = read_table6(token)
        checks.append(("L  table6.dat = the typed radii and mass fractions",
                       tab == list(zip(R_TAB, M_TAB))))
        nn = log.count("[NOTICE]")
        checks.append((f"W  one notice that column 6 lies outside the table ({nn} printed)",
                       nn == 1 and "the radius in column 6 lies outside" in log))
    fail = 0
    for name, ok in checks:
        print(f"[{'PASS' if ok else 'FAIL'}] {name}")
        fail += (not ok)
    print(f"       diameters: min {min(d):.4e} max {max(d):.4e} mean {sum(d) / n:.4e} m; "
          f"mass flow per parcel {w[0]:.6e} kg/s")
    return 1 if fail else 0


if __name__ == "__main__":
    sys.exit(main())
