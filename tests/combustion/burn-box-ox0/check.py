#!/usr/bin/env python3
"""burn-box-ox0: an oxide fraction at injection (xi-ox0) on the capped burning droplet, against an independent oracle.

The burning line carries xi-ox0 = x0 = 0.1. The inlet's diameter d_p0 and mass flow are the whole particle's, metal and
oxide at the equal volume: the parcel's mass is m_p0 = (pi/6) d_p0^3/((1 - x0)/rho_Al + x0/rho_res), its metal
m_Al0 = (1 - x0) m_p0, its cap m_ox0 = x0 m_p0 and its number rate mdot/m_p0. The cap grows from m_ox0 by xi nu of the
metal burned; the burn law acts on the metal's diameter, scaled by the exposed-surface factor F(f_ox) and by the
stretch re-evaluated from the initial cap,

  ds/dt = -K_eff c F(f_ox(s))/d_Al0^n,  s = (d_Al/d_Al0)^n,  m_Al = m_Al0 s^(3/n),  m_ox = m_ox0 + xi nu (m_Al0 - m_Al)
  c = (d_Al0/d_p0)^n int_0^1 ds/F(f_ox(s)) = 1.388383   (exact F, midpoint 20000 cells)

so the particle burns out in Beckstead's time for its injected diameter, t_b = d_p0^n/K_eff. At burnout the cap becomes
one parcel of the product OX (inert: its line carries no model). F comes from Gallier's lobe by bisection
(burn-box-cap's oracle), never from IGLOO's table; RK4 in time on (s, T) from (1, T_g) at t = 0 to each row's printed
t. Nothing is read from the code under test to build an expectation: the inputs below and the printed rows.

  W   the log: the xi-ox0 line, the stretch c read back within 3e-6 of 1.388383, the capped banner, 25 residues
  I   every burning parcel's injection row: m = m_p0 = 4.278788e-11 kg within 2e-6, d = d_p0
  C   every burning parcel's row before its exit record against the oracle at its t: m and d 4e-6, T 0.1 K; the
      coast u = u_g, v = w = 0
  B   each parent's exit record (bc 0) at a time inside [t(s_b), t_b] = [7.1974415e-03, 7.2246741e-03] s, 3e-6 at
      each end; s_b = (1e-15 kg/m_Al0)^(n/3), the burnout threshold on the metal
  R1  one residue per parent burned out inside the box: the product's IDs 1-25 = the parents' exit records with
      X < Lx - 1e-6; the product's banner counts them
  R2  each residue's birth row: m = m_ox0 + xi nu (m_Al0 - 1e-15 kg) = 1.883055e-11 kg within 3e-5, d = (6 m/(pi
      rho_res))^(1/3) of the printed m within 5e-6
  R3  sum of wdot(A) = sum of the parents' mdot - sum of the residues' mdot (exit records), wdot(OX) = 0: 1e-5
  R3e sum of the energy source = sum_parents mdot (h_mix(600 K) + u_g^2/2) - sum_residues mdot (c_res T_exit +
      hOff_res + u^2/2), h_mix = (1 - x0) h_A + x0 h_OX on each datum: 1e-5
  R5  each residue's mdot = (m_res/m_p0) x its parent's: 1e-5
  S   the scatter quanta of A and of OX = (sum_p mdot_p/m_p0) (L_ref/u_g)/(100 x 25), L_ref the box diagonal: 1e-3

X0 = 0 checks a run of the line without xi-ox0 (c 1.397686, m_p0 = (pi/6) rho_Al d_p0^3, no xi-ox0 line in the log).

The gate fails on the stretch not re-evaluated (W, C1 3.2e-3, B: the burnouts 0.63 % before the window), on the factor
(d_Al0/d_p0)^n dropped (W, C, B 4.6 % early), on the split on (pi/6) rho_Al d_p0^3 (I, R2, R5 2.5e-2), on the split
skipped (C1 0.14, R2 0.14, R3e 0.35), on the number rate left at mdot/((pi/6) rho_Al d_p0^3) (R3 4.7e-2, R5 2.5e-2), on
the scatter pre-pass's number rate on (pi/6) rho_Al d_p0^3 (S 2.6e-2) and on the xi-ox0 line missing from the log (W).
"""
import importlib.util
import math
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "tools"))
from datcols import zones  # noqa: E402

_spec = importlib.util.spec_from_file_location("burnboxcap_oracle", os.path.join(HERE, "..", "burn-box-cap", "check.py"))
cap = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(cap)                     # Gallier's lobe (exact F), burn-box's inputs (cap.bb), xi, nu
bb = cap.bb

# ---- inputs: phase.txt (both lines), the table's zones, bc.txt, the box ----
X0 = 0.1                                          # xi-ox0 on the burning line
STRETCH = {0.0: 1.397686, 0.1: 1.388383}          # (d_Al0/d_p0)^n int_0^1 ds/F(f_ox(s; x0)): exact F, 20000 cells
KEFF = 2.0e-6 * bb.XEFF
XI, NU = cap.XI, cap.NU
RHO_AL, C_AL, D0, NB = bb.RHO_P, bb.CP_P, bb.D0, bb.NB
RHO_RES, C_RES = 3950.0, 1887.6
HOFF_A, HOFF_RES = -372687.5, -1.723894e7
MB = 1.0e-15
LX, LY, LZ = 0.15, 0.05, 0.05
N_PARENTS, TRAJ_SAMPLE = 25, 100
DT = 1.0e-6

C_STRETCH = STRETCH[X0]
MP0 = math.pi / 6.0 * D0 ** 3 / ((1.0 - X0) / RHO_AL + X0 / RHO_RES)
MAL0, MOX0 = (1.0 - X0) * MP0, X0 * MP0
DAL0 = (6.0 * MAL0 / (math.pi * RHO_AL)) ** (1.0 / 3.0)
MRES = MOX0 + XI * NU * (MAL0 - MB)
SB = (MB / MAL0) ** (NB / 3.0)
TB = D0 ** NB / KEFF

TOL_C, TOL_I = 3.0e-6, 2.0e-6
TOL_M, TOL_D, TOL_T, TOL_U, TOL_VW = 4.0e-6, 4.0e-6, 0.1, 1.0e-4, 1.0e-6
TOL_B = 3.0e-6
TOL_M2, TOL_D2 = 3.0e-5, 5.0e-6
TOL_R3, TOL_R3E, TOL_R5, TOL_S = 1.0e-5, 1.0e-5, 1.0e-5, 1.0e-3
MIN_ROWS = 25


def state(s):
    """(m_Al, m_ox, d_Al, d) at s = (d_Al/d_Al0)^n, the cap from its initial mass."""
    s = max(s, 0.0)
    mal = MAL0 * s ** (3.0 / NB)
    mox = MOX0 + XI * NU * (MAL0 - mal)
    dal = DAL0 * s ** (1.0 / NB)
    deq = (6.0 / math.pi * (mal / RHO_AL + mox / RHO_RES)) ** (1.0 / 3.0)
    return mal, mox, dal, deq


def rhs(s, T):
    mal, mox, dal, deq = state(s)
    F = cap.LOBE.F(mox / (mal + mox))
    ds = -KEFF * C_STRETCH * F / DAL0 ** NB if s > 0.0 else 0.0
    mdot = RHO_AL * math.pi * KEFF / (2.0 * NB) * dal ** (3.0 - NB) * C_STRETCH * F
    dT = (2.0 * math.pi * deq * bb.KG * (bb.T_G - T) + bb.BETA * bb.QCOMB * mdot) / (mal * C_AL + mox * C_RES)
    return ds, dT


def rk4(s, T, h):
    k1 = rhs(s, T)
    k2 = rhs(s + 0.5 * h * k1[0], T + 0.5 * h * k1[1])
    k3 = rhs(s + 0.5 * h * k2[0], T + 0.5 * h * k2[1])
    k4 = rhs(s + h * k3[0], T + h * k3[1])
    return (s + h / 6.0 * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0]),
            T + h / 6.0 * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1]))


def oracle(times):
    """(s, T) at each of the sorted times and the time s crosses s_b (linear within its step): RK4 from (1, T_g) at
    t = 0 in steps of at most DT, continued past the last time until the crossing."""
    out, st = {}, {"s": 1.0, "T": bb.T_G, "t": 0.0, "t_sb": None}

    def step(h):
        s1, T1 = rk4(st["s"], st["T"], h)
        if st["t_sb"] is None and s1 <= SB < st["s"]:
            st["t_sb"] = st["t"] + h * (st["s"] - SB) / (st["s"] - s1)
        st["s"], st["T"], st["t"] = s1, T1, st["t"] + h

    for tq in times:
        while st["t"] < tq - 1e-15:
            step(min(DT, tq - st["t"]))
        out[tq] = (st["s"], st["T"])
    while st["t_sb"] is None and st["t"] < 2.0 * TB:
        step(DT)
    return out, st["t_sb"] if st["t_sb"] is not None else float("nan")


def read_source(path):
    """Cell sums of wdot(A), wdot(OX) and E from source.tec, and whether wdot(OX) is zero everywhere."""
    lines = open(path).read().splitlines()
    zi = next(i for i, l in enumerate(lines) if l.strip().upper().startswith("ZONE"))
    I, J, K = (int(re.search(r"%s=\s*(\d+)" % a, lines[zi]).group(1)) for a in "IJK")
    nn, nc = I * J * K, (I - 1) * (J - 1) * (K - 1)
    vals = [float(v) for l in lines[zi + 1:] for v in l.split()]
    if len(vals) != 3 * nn + 6 * nc:
        raise ValueError("%s: %d values, expected %d" % (path, len(vals), 3 * nn + 6 * nc))
    var = [vals[3 * nn + k * nc: 3 * nn + (k + 1) * nc] for k in range(6)]   # wdot(A) wdot(OX) Fx Fy Fz E
    return sum(var[0]), all(v == 0.0 for v in var[1]), sum(var[5])


def by_id(path):
    """Rows of a trajectory file per ID, sorted by t."""
    out = {}
    for z in zones(path):
        for r in z["rows"]:
            out.setdefault(r["ID"], []).append(r)
    for rs in out.values():
        rs.sort(key=lambda r: r["t"])
    return out


def log_checks(log):
    """W: the xi-ox0 line (none at X0 = 0), the stretch line read back, the capped banner, the product's banner."""
    banner = re.findall(r"combustion product \(residues: (\d+)\)", log)
    nres = int(banner[0]) if len(banner) == 1 else -1
    capped = "law with an oxide cap (xi-cap, product OX)" in log
    tok = re.findall(r"\] xi-ox0=\s*(\S+)", log)
    if X0 > 0.0:
        m = re.findall(r"\] oxide cap on product OX: stretch c = \s*([0-9.]+) from xi-ox0 = ([0-9.]+)\s*$", log, re.M)
        lines_ok = tok == ["%.3E" % X0] and len(m) == 1 and m[0][1] == "%6.4f" % X0
    else:
        m = [(v,) for v in re.findall(r"\] oxide cap on product OX: stretch c = \s*([0-9.]+)\s*$", log, re.M)]
        lines_ok = not tok and len(m) == 1
    c = float(m[0][0]) if len(m) == 1 else float("nan")
    err = abs(c / C_STRETCH - 1.0)
    return [("W the log: xi-ox0 line %s (want %s), stretch line %s; c %.6f vs %.6f (rel %.2e <= %.0e); the capped "
             "banner %s; the product's banner %d residues (need %d)"
             % (" ".join(tok) if tok else "none", "%.3E" % X0 if X0 > 0.0 else "none",
                "found" if len(m) == 1 else "MISSING", c, C_STRETCH, err, TOL_C, "found" if capped else "MISSING",
                nres, N_PARENTS),
             lines_ok and err <= TOL_C and capped and nres == N_PARENTS)], nres


def main():
    out = os.path.join(HERE, "OUTPUT")
    for f in ("trajectories-A.dat", "outloc-A.dat", "trajectories-OX.dat", "outloc-OX.dat", "source.tec"):
        if not os.path.exists(os.path.join(out, f)):
            print("[FAIL] OUTPUT/%s missing" % f)
            return 1
    log = open(os.path.join(HERE, "run_out.txt"), errors="replace").read()
    checks, nres = log_checks(log)

    rows = by_id(os.path.join(out, "trajectories-A.dat"))
    exits = [r for z in zones(os.path.join(out, "outloc-A.dat")) for r in z["rows"]]
    ex = {r["ID"]: r for r in exits}

    # I: the injection rows
    wi, bad_d, n_inj = 0.0, 0, 0
    for pid, rs in rows.items():
        r = rs[0]
        if r["t"] != 0.0:
            continue
        n_inj += 1
        wi = max(wi, abs(r["m<sub>p"] / MP0 - 1.0))
        bad_d += abs(r["d<sub>p"] / D0 - 1.0) > 1.0e-12
    checks.append(("I injection rows of %d parcels (need %d): m vs m_p0 = %.6e kg, worst %.2e <= %.0e; d = d_p0 %.1e "
                   "m: %d differ" % (n_inj, N_PARENTS, MP0, wi, TOL_I, D0, bad_d),
                   n_inj == N_PARENTS and wi <= TOL_I and bad_d == 0))

    # C: every row before its parcel's exit record against the oracle; B: the burnout times
    before = [(pid, r) for pid, rs in rows.items() for r in rs if pid in ex and r["t"] < ex[pid]["t"]]
    orc, t_sb = oracle(sorted({r["t"] for _, r in before}))
    wm = wd = wT = coast = vw = 0.0
    at_m = at_d = at_T = None
    for pid, r in before:
        s, T = orc[r["t"]]
        mal, mox, dal, deq = state(s)
        em = abs(r["m<sub>p"] - (mal + mox)) / (mal + mox)
        ed = abs(r["d<sub>p"] - deq) / deq
        eT = abs(r["T"] - T)
        if em > wm:
            wm, at_m = em, r["t"]
        if ed > wd:
            wd, at_d = ed, r["t"]
        if eT > wT:
            wT, at_T = eT, r["t"]
        coast = max(coast, abs(r["U"] - bb.U_G))
        vw = max(vw, abs(r["V"]), abs(r["W"]))
    nmin = min((sum(1 for p, _ in before if p == pid) for pid in ex), default=0)
    checks += [
        ("C0 non-vacuity: %d rows before the exit records of %d parcels, at least %d per parcel (need %d parcels, %d "
         "rows each)" % (len(before), len(ex), nmin, N_PARENTS, MIN_ROWS), len(ex) == N_PARENTS and nmin >= MIN_ROWS),
        ("C1 m = m_Al + m_ox vs the oracle: worst %.3e <= %.0e%s"
         % (wm, TOL_M, " at t = %.6e s" % at_m if at_m is not None else ""), wm <= TOL_M),
        ("C2 d (equal-volume sphere) vs the oracle: worst %.3e <= %.0e%s"
         % (wd, TOL_D, " at t = %.6e s" % at_d if at_d is not None else ""), wd <= TOL_D),
        ("C3 T vs the oracle: worst %.3e K <= %.1f K%s"
         % (wT, TOL_T, " at t = %.6e s" % at_T if at_T is not None else ""), wT <= TOL_T),
        ("C4 the coast: |u - u_g| %.2e <= %.0e, |v|,|w| %.2e <= %.0e" % (coast, TOL_U, vw, TOL_VW),
         coast <= TOL_U and vw <= TOL_VW),
    ]
    tex = sorted(r["t"] for r in exits)
    lo, hi = t_sb * (1.0 - TOL_B), TB * (1.0 + TOL_B)
    nb0 = sum(1 for r in exits if r["bc"] == 0)
    checks.append(("B burnout of every parent: %d exit records, %d with bc 0; t in [%.7e, %.7e] s against the window "
                   "[t(s_b), t_b] = [%.7e, %.7e] s (%.0e at each end)"
                   % (len(exits), nb0, tex[0] if tex else 0.0, tex[-1] if tex else 0.0, t_sb, TB, TOL_B),
                   len(exits) == N_PARENTS and nb0 == N_PARENTS and tex[0] >= lo and tex[-1] <= hi))

    # R1, R2: the residues
    res = by_id(os.path.join(out, "trajectories-OX.dat"))
    res_exit = {r["ID"]: r for z in zones(os.path.join(out, "outloc-OX.dat")) for r in z["rows"]}
    inside = [r for r in exits if r["X"] < LX - 1e-6]
    ids = sorted(res)
    checks.append(("R1 one residue per parent burned out inside the box: %d residues (IDs %s), %d parents inside, "
                   "%d exit records, banner %d (need %d)"
                   % (len(ids), "1-%d" % max(ids) if ids else "none", len(inside), len(res_exit), nres, N_PARENTS),
                   len(ids) == len(inside) == len(res_exit) == nres == N_PARENTS
                   and ids == list(range(1, N_PARENTS + 1)) and sorted(res_exit) == ids))
    wm2 = wd2 = 0.0
    for rs in res.values():
        m, d = rs[0]["m<sub>p"], rs[0]["d<sub>p"]
        wm2 = max(wm2, abs(m - MRES) / MRES)
        wd2 = max(wd2, abs(d - (6.0 * m / (math.pi * RHO_RES)) ** (1.0 / 3.0)) / d)
    checks.append(("R2 residue mass vs m_ox0 + xi nu (m_Al0 - 1e-15 kg) = %.6e kg: worst %.2e <= %.0e; d vs "
                   "(6 m/(pi rho_res))^(1/3) of the printed m: worst %.2e <= %.0e" % (MRES, wm2, TOL_M2, wd2, TOL_D2),
                   len(res) > 0 and wm2 <= TOL_M2 and wd2 <= TOL_D2))

    # R3, R3e: the books
    wA, ox_zero, E = read_source(os.path.join(out, "source.tec"))
    mpar = sum(r["mdot"] for r in exits)
    mres = sum(r["mdot"] for r in res_exit.values())
    exp_m = mpar - mres
    r3 = abs(wA - exp_m) / exp_m
    h_mix = (1.0 - X0) * (C_AL * bb.T_G + HOFF_A) + X0 * (C_RES * bb.T_G + HOFF_RES)
    e_in = mpar * (h_mix + 0.5 * bb.U_G ** 2)
    e_out = sum(r["mdot"] * (C_RES * r["T"] + HOFF_RES + 0.5 * r["|u<sub>p</sub>|"] ** 2) for r in res_exit.values())
    exp_e = e_in - e_out
    r3e = abs(E - exp_e) / abs(exp_e)
    checks += [
        ("R3 sum wdot(A) = sum mdot(parents) - sum mdot(residues): %.9e vs %.9e kg/s, rel %.2e <= %.0e; wdot(OX) = 0"
         % (wA, exp_m, r3, TOL_R3), r3 <= TOL_R3 and ox_zero),
        ("R3e sum E = parents' inflow at h_mix(600 K) %.6e J/kg - residues' outflow: %.6e vs %.6e W, rel %.2e <= %.0e"
         % (h_mix, E, exp_e, r3e, TOL_R3E), r3e <= TOL_R3E),
    ]

    # R5: the residue's flow
    w5 = 0.0
    for rid, r in res_exit.items():
        if rid not in res or rid not in ex:
            w5 = float("inf")
            continue
        w5 = max(w5, abs(r["mdot"] - res[rid][0]["m<sub>p"] / MP0 * ex[rid]["mdot"]) / r["mdot"])
    checks.append(("R5 residue mdot = (m_res/m_p0) x parent mdot: worst %.2e <= %.0e" % (w5, TOL_R5),
                   len(res_exit) > 0 and w5 <= TOL_R5))

    # S: the scatter quanta
    lref = math.sqrt(LX ** 2 + LY ** 2 + LZ ** 2)
    exp_s = sum(r["mdot"] / MP0 for r in exits) * (lref / bb.U_G) / (TRAJ_SAMPLE * N_PARENTS)
    got = {}
    for mat in ("A", "OX"):
        v = re.findall(r"scatter cloud \[%s\]: dNscat=\s*(\S+)" % mat, log)
        got[mat] = float(v[0]) if len(v) == 1 else float("nan")
    ws = max(abs(got[k] / exp_s - 1.0) for k in got)
    checks.append(("S scatter quanta dNscat [A] %.3e, [OX] %.3e vs (sum mdot_p/m_p0)(L_ref/u_g)/(%d x %d) = %.4e: "
                   "worst rel %.2e <= %.0e" % (got["A"], got["OX"], TRAJ_SAMPLE, N_PARENTS, exp_s, ws, TOL_S),
                   ws <= TOL_S))

    print("oracle: c %.6f, m_p0 %.6e kg, m_Al0 %.6e kg, m_ox0 %.6e kg, d_Al0 %.6e m, m_res %.6e kg; t(s_b) %.7e s, "
          "t_b %.7e s" % (C_STRETCH, MP0, MAL0, MOX0, DAL0, MRES, t_sb, TB))
    fail = 0
    for name, ok in checks:
        print("[%s] %s" % ("PASS" if ok else "FAIL", name))
        fail += (not ok)
    print("\n[%s] burn-box-ox0: the oxide fraction at injection %s the independent oracle"
          % ("PASS" if fail == 0 else "FAIL", "matches" if fail == 0 else "misses"))
    return 1 if fail else 0


if __name__ == "__main__":
    sys.exit(main())
