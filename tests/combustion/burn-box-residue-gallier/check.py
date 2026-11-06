#!/usr/bin/env python3
"""burn-box-residue-gallier: the burnout residue under Gallier et al. 2021's deposition law, against the law's path
and against Gallier's own residue sizes.

Under the law the cap's growth in the metal mass, dm_ox/dm_Al = -nu xi(d_Al, f_ox), d_Al = d_Al0 (m_Al/m_Al0)^(1/3),
holds no time, burn law or stretch: the residue's mass over the injected one, z = m_ox/m0 at burnout, depends on the
law and the injected size alone (law_path.path, the law built from the committed tables, never from IGLOO). Five
inlet rows carry 20, 70, 200, 400 and 100 um (burn-box's bc.txt with the radius set by the row j); K-burn 1.5e-4
burns every parcel out inside the box. A parcel burns out at the first accepted step with m_Al <= m_b,
m_b = min(1e-6 m0, 1e-20 kg) (tools/burnout.py: 1e-20 kg at every size here), so its residue lies between the path's
cap at m_Al = m_b and at m_Al = 0, and its burnout time between the path's time to m_Al = m_b and Beckstead's
d0^n/K_eff (the stretch per parcel keeps every burn time on it).

  R1  one residue per parcel, every parcel burned out inside the box (exit X < Lx - 1e-6); the banner counts 25
  R2  each residue's mass at birth (its first row's m) over m0 = rho_Al pi d0^3/6, d0 the inlet size (the parent's
      first row's d to the micron): inside [z(m_b), z(0)] of the path from d0, each end widened by 1e-5 relative
      (the E13.6 print of m, the integrator, the last step's overshoot past m_Al = 0: its metal clamped, its cap kept)
  T1  each parent's burnout time (its exit record's t): inside [t(m_b), d0^n/K_eff], t(m_b) the path's fraction
      of d0^n/K_eff to m_Al = m_b, each end widened by 1e-5 relative (the printed t, the integrator)
  R3  sum of wdot(A) = sum of the parents' mdot - sum of the residues' mdot, wdot(OX) = 0; R3e the energy books
      (burn-box-residue's books): 1e-5
  B   Gallier's residue sizes: beta = (z rho_Al/rho_ox)^(1/3) with Gallier's densities at 2570 K (1673.7 and 2820.3:
      z does not depend on the densities, beta does), from the measured residue mass of each size, against
      Gallier's 0.57 (20 um), 0.55 (70 um, Table 3 at f_ox0 = 0), 0.49 (200 um) and 0.40 (400 um, PDF p. 14):
      0.03 (their two digits and the digitization)

The gate fails on the constant law at h(0) (xi-cap 0.069522: R2 25 of 25, B beta 0.427 at every size), on g read at
the injected diameter (R2 and T1 20 of 25, B beta 0.256 at 400 um), on the law read at the drop's diameter (R2 and
T1 20 of 25) and on the stretch left at 1 (T1 25 of 25: every burnout 1.20-1.34 times Beckstead's time). The 400 um
size is the one B needs against g at the injected diameter (at 200 um that mutant's beta, 0.460, passes B).
"""
import importlib.util
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "tools"))
from datcols import zones  # noqa: E402
from burnout import burn_mass  # noqa: E402


def _load(name, rel):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, rel))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


law = _load("law_path", os.path.join("..", "burn-box-cap-gallier", "law_path.py"))
res = _load("burnboxresidue", os.path.join("..", "burn-box-residue", "check.py"))   # the books (R3, R3e)
bb = law.cap.bb

RHO_AL, RHO_RES, NB = bb.RHO_P, 3950.0, bb.NB
KEFF = 1.5e-4 * bb.XEFF
LX = 0.15
N_PARENTS = 25
RHO_AL_G, RHO_OX_G = 3226.0 - 0.604 * 2570.0, 3056.0 - 0.97 * (2570.0 - 2327.0)
GALLIER = {20: 0.57, 70: 0.55, 200: 0.49, 400: 0.40}
TOL_M, TOL_T, TOL_B = 1.0e-5, 1.0e-5, 0.03


def main():
    if not os.path.exists("OUTPUT/outloc-OX.dat"):
        print("[FAIL] OUTPUT/outloc-OX.dat missing")
        return 1
    par, res_exit, res_rows, wA, ox_zero, E, log = res.load_run(".")
    parents = [r for g in par for r in g]
    first = {}
    for z in zones("OUTPUT/trajectories-A.dat"):
        for r in z["rows"]:
            first.setdefault(r["ID"], r)
    inside = [p for p in parents if p["X"] < LX - 1e-6]
    checks = [
        ("R1 one residue per parcel burned out inside the box: %d parents, %d inside, residue IDs %s"
         % (len(parents), len(inside), "= theirs" if sorted(res_rows) == sorted(p["ID"] for p in inside) else "differ"),
         len(parents) == N_PARENTS and len(inside) == N_PARENTS and sorted(res_rows) == sorted(p["ID"] for p in inside)
         and "combustion product (residues: %d)" % N_PARENTS in log),
    ]

    paths, bad_m, bad_t, worst = {}, [], [], {}
    for p in sorted(parents, key=lambda q: q["ID"]):
        key = round(first[p["ID"]]["d<sub>p"] * 1e6)          # the inlet size, its first row's d to the micron
        d0 = key * 1.0e-6
        m0 = RHO_AL * math.pi / 6.0 * d0 ** 3
        if key not in paths:
            c, z0 = law.path(d0, NB, 0.0, RHO_AL, RHO_RES)
            _, zb, fb = law.path(d0, NB, 0.0, RHO_AL, RHO_RES, x_end=burn_mass(m0) / m0, frac=True)
            paths[key] = (z0, zb, fb, d0 ** NB / KEFF)
        z0, zb, fb, tb = paths[key]
        rows = res_rows.get(p["ID"], [])
        if not rows:
            bad_m.append(p["ID"])
            continue
        z = rows[0]["m<sub>p"] / m0
        if not zb * (1.0 - TOL_M) <= z <= z0 * (1.0 + TOL_M):
            bad_m.append(p["ID"])
        if not tb * fb * (1.0 - TOL_T) <= p["t"] <= tb * (1.0 + TOL_T):
            bad_t.append(p["ID"])
        w = worst.setdefault(key, [z, z, p["t"], p["t"]])
        w[0], w[1], w[2], w[3] = min(w[0], z), max(w[1], z), min(w[2], p["t"]), max(w[3], p["t"])
    for key in sorted(paths):
        z0, zb, fb, tb = paths[key]
        zl, zh, tl, th = worst[key]
        print("d0 %3d um: z in [%.9f, %.9f] (path %.9f at m_b, %.9f at 0); t in [%.6e, %.6e] s (path %.6e at m_b, "
              "d0^n/K_eff %.6e)" % (key, zl, zh, zb, z0, tl, th, tb * fb, tb))
    checks += [
        ("R2 each residue's mass over m0 inside [z(m_b), z(0)] of the law's path (+-%.0e): %d of %d outside"
         % (TOL_M, len(bad_m), len(parents)), not bad_m and len(paths) == 5),
        ("T1 each burnout time inside [t(m_b), d0^n/K_eff] (+-%.0e): %d of %d outside" % (TOL_T, len(bad_t),
                                                                                       len(parents)), not bad_t),
    ]
    checks += res.books("", par, res_exit, wA, ox_zero, E)

    bres = []
    for key, gv in sorted(GALLIER.items()):
        if key not in worst:
            bres.append((key, gv, float("nan")))
            continue
        z = 0.5 * (worst[key][0] + worst[key][1])
        bres.append((key, gv, (z * RHO_AL_G / RHO_OX_G) ** (1.0 / 3.0)))
    worst_b = max(abs(b - gv) for _, gv, b in bres)
    checks.append(("B Gallier's residue sizes beta = (z rho_Al/rho_ox)^(1/3) at 2570 K: %s; worst %.3f <= %.2f" % (
        ", ".join("%d um %.3f (%.2f)" % (k, b, gv) for k, gv, b in bres), worst_b, TOL_B), worst_b <= TOL_B))

    fail = 0
    for name, ok in checks:
        print("[%s] %s" % ("PASS" if ok else "FAIL", name.strip()))
        fail += (not ok)
    print("\n[%s] burn-box-residue-gallier: the residues under Gallier's law %s the law's path and Gallier's sizes"
          % ("PASS" if fail == 0 else "FAIL", "match" if fail == 0 else "miss"))
    return 1 if fail else 0


if __name__ == "__main__":
    sys.exit(main())
