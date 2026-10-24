#!/usr/bin/env python3
"""burn-box-residue: the burnout residue of an aluminium droplet against independent books.

A capped parcel (ODE model 8, xi-cap 0.2, residue = on) burns out near x = 0.072 m. Its cap becomes one parcel of the
combustion product OX (a 0-group line of the phase file, its own table zone and solidification data), born at the
parent's position, velocity, temperature and time with the cap's mass, written as a birth row; the gas books lose the
residue's birth flux at that point. OX freezes: born undercooled (f = 0), it nucleates at its first segment (the
whole drop freezes, enthalpy kept) and relaxes to the gas temperature. Nothing is read from the code under test to
build an expectation: the inputs below, the product's line and table, and the printed rows.

  R1  one residue per parcel that burned out inside the box: the product's IDs = the parents' exit records with
      X < Lx - 1e-6 (25 here); the product's banner counts them
  R2  each residue's birth row = its parent's exit record (X, Y, Z, T printed alike, the speed within its print);
      m = xi nu (m0 - 1e-15 kg) within 3e-5 (the burnout threshold is 2.4e-5 of m0), d = (6 m/(pi rho_res))^(1/3)
      of the printed m within 5e-6 (the E13.6 prints of d and m)
  R3  sum of wdot(A) = sum of the parents' mdot - sum of the residues' mdot (exit records), wdot(OX) = 0: 1e-5
  R3e sum of the energy source = sum_parents mdot (h_A(600 K) + u^2/2) - sum_residues mdot (h_OX(T_exit) + u^2/2),
      h_OX on the product's datum by phase (f of the exit record): 1e-5
  R4  each residue's rows: born undercooled with f = 0, solid (f = 1) after its first segment at the temperature of
      the enthalpy-conserving jump, then T = T_g + (T_a - T_g) exp(-(t - t_b)/tau), tau = rho_res c_s d^2/(12 k_g):
      the print budget of temp-relax with the jump's anchor (c_l/c_s) and d's print in tau
  R5  each residue's mdot = (m_res/m0) x its parent's: 1e-5
  D   ref-dres/: d-residue = 5e-6 m and the burning line in two groups (bc.txt's face table written twice): 50
      distinct residue IDs, 1-25 and 26-50 (the second group's at 25 + its parents' IDs); d = 5e-6 m; npdot (the exit
      record's mdot over the birth row's m) = 60.958139 x its parent's (4e-5); R3, R3e; T = 600 K at the outlet
  F   ref-fail/: [IGLOO-ODE] fail-after = 3.6e-3 s: every parent ends at t* >= 3.6e-3 s; its exit record is its
      residue's birth row; m_res = xi nu (m0 - m_Al(t*)), m_Al(t*) from the burn law (Gallier's F by bisection,
      c = 1.397686) at the printed t*: 1e-5; R1, R3, R3e

The gate fails on a residue sized with the burning material's density (R2: d 2.168897e-5 m), on a birth flux left
in the gas books (R3: 0.378 of the injected flow), on the product's datum left out of the birth flux (R3e 100.7 %),
on a residue born solid at its parent's temperature (R4: 63.9 K at the outlet), on d-residue ignored (D), on a
residue ID without its group's offset (D: 25 distinct IDs) and on a failure exit that keeps the segment-start
state (F: the residue's mass).
"""
import importlib.util
import math
import os
import re
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "tools"))
from datcols import zones  # noqa: E402

_spec = importlib.util.spec_from_file_location("burnboxcap_oracle", os.path.join(HERE, "..", "burn-box-cap", "check.py"))
cap = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(cap)                     # Gallier's lobe (exact F), burn-box's inputs (cap.bb), xi, nu, c
bb = cap.bb
IGLOO = os.environ.get("IGLOO_EXE", os.path.abspath(os.path.join(HERE, "..", "..", "..", "bin", "IGLOO")))

# ---- inputs: phase.txt (both lines) and the table's zones ----
KEFF = 2.0e-6 * bb.XEFF
XI, NU, C_STRETCH = cap.XI, cap.NU, cap.C_STRETCH
RHO_RES, CL, CS, HF, TM = 3950.0, 1887.6, 1360.0, 1.08949e6, 2327.0
TN = 0.8 * TM
HOFF_A, HOFF_L = -372687.5, -1.723894e7
M0 = cap.M0
MB = 1.0e-15
LX = 0.15
D_RES = 5.0e-6
FAIL_AFTER = 3.6e-3
N_PARENTS = 25

TOL_M2, TOL_D2, TOL_V = 3.0e-5, 5.0e-6, 1.5e-6
TOL_R3, TOL_R3E, TOL_R5, TOL_DNP, TOL_F = 1.0e-5, 1.0e-5, 1.0e-5, 4.0e-5, 1.0e-5
HALF_ULP, INT_FLOOR, PRINT_D = 0.5e-6, 1.0e-8, 5.0e-7


def h_ox(T, f):
    """The product's specific enthalpy on its datum by phase: liquid or undercooled, plateau, solid."""
    if f >= 1.0:
        return CL * TM + HOFF_L - HF - CS * (TM - T)
    if f > 0.0:
        return CL * TM + HOFF_L - f * HF
    return CL * T + HOFF_L


def jump(Tb):
    """Nucleation of an undercooled drop at Tb, the enthalpy kept: (T, f) after it."""
    f0 = CL * (TM - Tb) / HF
    if f0 < 1.0:
        return TM, f0
    return TM - (CL * (TM - Tb) - HF) / CS, 1.0


def m_al_at(t):
    """The metal's mass at time t from the burn law with the cap: ds/dt = -K_eff c F(f_ox(s))/d0^n, RK4."""
    def rate(s):
        x = max(s, 0.0) ** (3.0 / bb.NB)
        f = XI * NU * (1.0 - x) / (x + XI * NU * (1.0 - x))
        return -KEFF * C_STRETCH * cap.LOBE.F(f) / bb.D0 ** bb.NB if s > 0.0 else 0.0
    s, tt, dt = 1.0, 0.0, 1.0e-6
    while tt < t - 1e-15:
        h = min(dt, t - tt)
        k1 = rate(s); k2 = rate(s + 0.5 * h * k1); k3 = rate(s + 0.5 * h * k2); k4 = rate(s + h * k3)
        s += h / 6.0 * (k1 + 2 * k2 + 2 * k3 + k4)
        tt += h
    return M0 * max(s, 0.0) ** (3.0 / bb.NB)


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


def load_run(d):
    """The run's parent exit records (per group), the residues' rows (per ID) and exit records, its source sums."""
    out = os.path.join(d, "OUTPUT")
    par = [z["rows"] for z in zones(os.path.join(out, "outloc-A.dat"))]
    res_exit = {r["ID"]: r for z in zones(os.path.join(out, "outloc-OX.dat")) for r in z["rows"]}
    res_rows = {}
    for z in zones(os.path.join(out, "trajectories-OX.dat")):
        for r in z["rows"]:
            res_rows.setdefault(r["ID"], []).append(r)
    for rs in res_rows.values():
        rs.sort(key=lambda r: r["t"])
    wA, ox_zero, E = read_source(os.path.join(out, "source.tec"))
    log = open(os.path.join(d, "run_out.txt"), errors="replace").read()
    return par, res_exit, res_rows, wA, ox_zero, E, log


def books(name, par, res_exit, wA, ox_zero, E):
    """R3 and R3e over every group of the run."""
    mpar = sum(r["mdot"] for g in par for r in g)
    mres = sum(r["mdot"] for r in res_exit.values())
    exp_m = mpar - mres
    r3 = abs(wA - exp_m) / exp_m
    e_in = mpar * (bb.CP_P * bb.T_G + HOFF_A + 0.5 * bb.U_G ** 2)
    e_out = sum(r["mdot"] * (h_ox(r["T"], r["f"]) + 0.5 * r["|u<sub>p</sub>|"] ** 2) for r in res_exit.values())
    exp_e = e_in - e_out
    r3e = abs(E - exp_e) / abs(exp_e)
    return [
        ("%s R3 sum wdot(A) = sum mdot(parents) - sum mdot(residues): %.9e vs %.9e kg/s, rel %.2e <= %.0e; wdot(OX) = 0"
         % (name, wA, exp_m, r3, TOL_R3), r3 <= TOL_R3 and ox_zero),
        ("%s R3e sum E = parents' inflow - residues' outflow at the product's datum: %.6e vs %.6e W "
         "(per kg injected %.6e J/kg), rel %.2e <= %.0e" % (name, E, exp_e, exp_e / mpar, r3e, TOL_R3E),
         r3e <= TOL_R3E),
    ]


def birth_identity(name, par_by_id, res_rows):
    """Each residue's first row against its parent's exit record: X, Y, Z, T alike, the speed within its print."""
    worst_v, bad = 0.0, 0
    for rid, rs in res_rows.items():
        p = par_by_id(rid)
        b = rs[0]
        if p is None or any(b[k] != p[k] for k in ("X", "Y", "Z", "T")):
            bad += 1
            continue
        worst_v = max(worst_v, abs(math.sqrt(b["U"] ** 2 + b["V"] ** 2 + b["W"] ** 2) - p["|u<sub>p</sub>|"]))
    return ("%s R2 every birth row = its parent's exit record (X, Y, Z, T alike; speed within %.1e: %.1e): %d of %d "
            "differ" % (name, TOL_V, worst_v, bad, len(res_rows)), bad == 0 and worst_v <= TOL_V and len(res_rows) > 0)


def banner(log):
    m = re.findall(r"combustion product \(residues: (\d+)\)", log)
    return int(m[0]) if len(m) == 1 else -1


def check_main(checks):
    par, res_exit, res_rows, wA, ox_zero, E, log = load_run(HERE)
    exits = par[0]
    inside = [r for r in exits if r["X"] < LX - 1e-6]
    ids = sorted(res_rows)
    checks.append(("R1 one residue per parent burned out inside the box: %d residues (IDs %s), %d parents inside, "
                   "banner %d (need %d)" % (len(ids), "1-%d" % max(ids) if ids else "none", len(inside), banner(log),
                                            N_PARENTS),
                   len(ids) == len(inside) == N_PARENTS == len(res_exit) == banner(log)
                   and ids == list(range(1, N_PARENTS + 1))))
    by_id = {r["ID"]: r for r in exits}
    checks.append(birth_identity("main", lambda rid: by_id.get(rid), res_rows))
    mres = XI * NU * (M0 - MB)
    wm = wd = 0.0
    for rs in res_rows.values():
        m, d = rs[0]["m<sub>p"], rs[0]["d<sub>p"]
        wm = max(wm, abs(m - mres) / mres)
        wd = max(wd, abs(d - (6.0 * m / (math.pi * RHO_RES)) ** (1.0 / 3.0)) / d)
    checks.append(("R2 residue mass vs xi nu (m0 - 1e-15 kg) = %.6e kg: worst %.2e <= %.0e; d vs (6 m/(pi rho_res))^(1/3)"
                   " of the printed m: worst %.2e <= %.0e" % (mres, wm, TOL_M2, wd, TOL_D2),
                   wm <= TOL_M2 and wd <= TOL_D2))
    checks += books("main", par, res_exit, wA, ox_zero, E)
    # R4: the residue's temperature history
    worst, nrow, phase_ok = 0.0, 0, True
    for rs in res_rows.values():
        b = rs[0]
        phase_ok = phase_ok and b["f"] == 0.0 and b["T"] < TN
        Ta, fa = jump(b["T"])
        tau = RHO_RES * CS * b["d<sub>p"] ** 2 / (12.0 * bb.KG)
        L = bb.U_G * tau
        for r in rs[1:]:
            phase_ok = phase_ok and r["f"] == fa == 1.0
            dt = r["t"] - b["t"]
            Tp = bb.T_G + (Ta - bb.T_G) * math.exp(-dt / tau)
            tol = (HALF_ULP * (1.0 + CL / CS + 2.0 * abs(r["T"] - bb.T_G) / L)
                   + 2.0 * PRINT_D * abs(r["T"] - bb.T_G) * dt / tau + INT_FLOOR)
            worst = max(worst, abs(r["T"] - Tp) / tol)
            nrow += 1
    Tex = sorted(r["T"] for r in res_exit.values())
    checks.append(("R4 residue T(x): born undercooled (f = 0), solid after the jump, then the Nu = 2 relaxation over "
                   "%d rows: worst resid/tol %.3f <= 1; exit T %.4f K" % (nrow, worst, Tex[0] if Tex else 0.0),
                   phase_ok and nrow >= 25 * 20 and worst <= 1.0))
    # R5: the residue's flow
    w5 = 0.0
    for rid, r in res_exit.items():
        m = res_rows[rid][0]["m<sub>p"]
        w5 = max(w5, abs(r["mdot"] - m / M0 * by_id[rid]["mdot"]) / r["mdot"])
    checks.append(("R5 residue mdot = (m_res/m0) x parent mdot: worst %.2e <= %.0e" % (w5, TOL_R5), w5 <= TOL_R5))


def sibling(name, phase_text, bc_copies, ode_extra):
    """A sibling run of this case's input into HERE/<name>/ with its own phase file, bc.txt copies, [IGLOO-ODE] lines."""
    d = os.path.join(HERE, name)
    shutil.rmtree(d, ignore_errors=True)
    os.makedirs(os.path.join(d, "OUTPUT"))
    os.makedirs(os.path.join(d, "INPUT"))
    for f in ("solfile.tec", "properties.dat"):
        os.symlink(os.path.realpath(os.path.join(HERE, "INPUT", f)), os.path.join(d, "INPUT", f))
    bc = open(os.path.join(HERE, "INPUT", "bc.txt")).read()
    open(os.path.join(d, "INPUT", "bc.txt"), "w").write(bc * bc_copies)
    open(os.path.join(d, "INPUT", "phase.txt"), "w").write(phase_text)
    ini = open(os.path.join(HERE, "input.ini")).read()
    if ode_extra:
        ini = ini.replace("[IGLOO-ODE]\n", "[IGLOO-ODE]\n" + ode_extra + "\n", 1)
    open(os.path.join(d, "input.ini"), "w").write(ini)
    with open(os.path.join(d, "run_out.txt"), "w") as out, open(os.path.join(d, "run_err.txt"), "w") as err:
        rc = subprocess.run([IGLOO], cwd=d, stdout=out, stderr=err).returncode
    return d, rc


def check_d(checks):
    phase = open(os.path.join(HERE, "INPUT", "phase.txt")).read()
    phase = phase.replace("\nA 1 ", "\nA 2 ", 1).replace("residue=on", "residue=on d-residue=5.0e-6", 1)
    d, rc = sibling("ref-dres", phase, 2, None)
    if rc != 0:
        checks.append(("D sibling run (ref-dres) exit %d" % rc, False))
        return
    par, res_exit, res_rows, wA, ox_zero, E, log = load_run(d)
    ids = sorted(res_rows)
    g1, g2 = ({r["ID"]: r for r in par[0]}, {r["ID"]: r for r in par[1]}) if len(par) == 2 else ({}, {})
    checks.append(("D two groups of %s parents; %d residue exit records, %d distinct IDs, banner %d (need 50: 1-25 "
                   "and 26-50)" % ("/".join(str(len(g)) for g in par), len(res_exit), len(set(ids)), banner(log)),
                   len(par) == 2 and len(g1) == len(g2) == N_PARENTS and ids == list(range(1, 2 * N_PARENTS + 1))
                   and len(res_exit) == 2 * N_PARENTS and banner(log) == 2 * N_PARENTS))
    checks.append(birth_identity("D", lambda rid: g1.get(rid) if rid <= N_PARENTS else g2.get(rid - N_PARENTS),
                                 res_rows))
    dd = max(abs(r["d<sub>p"] - D_RES) for rs in res_rows.values() for r in rs) if res_rows else 1.0
    ratio_ref = XI * NU * (M0 - MB) / (RHO_RES * math.pi / 6.0 * D_RES ** 3)
    wn = 0.0
    for rid, r in res_exit.items():
        p = g1[rid] if rid <= N_PARENTS else g2[rid - N_PARENTS]
        ratio = (r["mdot"] / res_rows[rid][0]["m<sub>p"]) / (p["mdot"] / M0)
        wn = max(wn, abs(ratio - ratio_ref) / ratio_ref)
    Tex = [r["T"] for r in res_exit.values()]
    checks.append(("D every residue at d-residue %.1e m (worst |d - d_res| %.1e); npdot x %.6f its parent's: worst "
                   "%.2e <= %.0e; exit T within 1e-4 K of %.1f K (%.6f..%.6f)"
                   % (D_RES, dd, ratio_ref, wn, TOL_DNP, bb.T_G, min(Tex), max(Tex)),
                   dd <= 0.5e-11 and wn <= TOL_DNP and all(abs(T - bb.T_G) <= 1e-4 for T in Tex)))
    checks += books("D", par, res_exit, wA, ox_zero, E)


def check_f(checks):
    d, rc = sibling("ref-fail", open(os.path.join(HERE, "INPUT", "phase.txt")).read(), 1,
                    "fail-after   = %.1e" % FAIL_AFTER)
    if rc != 0:
        checks.append(("F sibling run (ref-fail) exit %d" % rc, False))
        return
    par, res_exit, res_rows, wA, ox_zero, E, log = load_run(d)
    exits = par[0]
    by_id = {r["ID"]: r for r in exits}
    tstar = [r["t"] for r in exits]
    checks.append(("F every parent ends at t* >= fail-after %.1e s: t* in [%.8e, %.8e]; %d residues, banner %d"
                   % (FAIL_AFTER, min(tstar), max(tstar), len(res_rows), banner(log)),
                   len(exits) == N_PARENTS and min(tstar) >= FAIL_AFTER and len(res_rows) == N_PARENTS
                   and banner(log) == N_PARENTS and sorted(res_rows) == list(range(1, N_PARENTS + 1))))
    checks.append(birth_identity("F", lambda rid: by_id.get(rid), res_rows))
    if not res_rows:
        checks.append(("F residue mass vs xi nu (m0 - m_Al(t*)): no residue born after the failure", False))
        return
    cache, wm = {}, 0.0
    for rid, rs in res_rows.items():
        t = by_id[rid]["t"]
        if t not in cache:
            cache[t] = XI * NU * (M0 - m_al_at(t))
        wm = max(wm, abs(rs[0]["m<sub>p"] - cache[t]) / cache[t])
    t0 = sorted(cache)[0]
    checks.append(("F residue mass vs xi nu (m0 - m_Al(t*)) (m_Al(t*)/m0 %.4f at t* = %.6e s): worst %.2e <= %.0e"
                   % (m_al_at(t0) / M0, t0, wm, TOL_F), wm <= TOL_F))
    checks += books("F", par, res_exit, wA, ox_zero, E)


def main():
    for f in ("outloc-A.dat", "outloc-OX.dat", "trajectories-OX.dat", "source.tec"):
        if not os.path.exists(os.path.join(HERE, "OUTPUT", f)):
            print("[FAIL] OUTPUT/%s missing" % f)
            return 1
    checks = []
    check_main(checks)
    check_d(checks)
    check_f(checks)
    fail = 0
    for name, ok in checks:
        print("[%s] %s" % ("PASS" if ok else "FAIL", name))
        fail += (not ok)
    print("\n[%s] burn-box-residue: the burnout residue %s its books"
          % ("PASS" if fail == 0 else "FAIL", "closes" if fail == 0 else "misses"))
    return 1 if fail else 0


if __name__ == "__main__":
    sys.exit(main())
