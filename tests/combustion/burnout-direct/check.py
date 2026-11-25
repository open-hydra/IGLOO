#!/usr/bin/env python3
"""burnout-direct: drops of combustion=direct heated into their burning temperature burn on the d^2-law and are removed
at d-burnout, inside the box, their remnant handed to the gas.

Inputs (not read from the solver's output): burn-box's uniform gas, T_g = 600 K, u_g = 10 m/s, k_g = 0.026 W/m/K,
gamma = 1.4, R = 287 J/kg/K, so c_pg = gamma R/(gamma - 1) = 1004.5 J/kg/K; 25 drops, d0 = 10 to 29.2 um, injected at
T0 = 300 K with the gas velocity (bc.txt kv = 1, kt = 0.5); rho_p = 2950 kg/m^3; the models row T-burn = 302 K,
q-comb = 3.1e7 J/kg, d-burnout = 3 um; Lv = 1.08e7 J/kg ([IGLOO-Properties]). Two runs: this directory, with
common/properties.dat (c_p 1250 J/kg/K constant: the temperature is the state), and the sibling ref-varcp/, the same
inputs with standard/temp-relax-varcp's table (c_p = 1000 + (T - 250) J/kg/K, rho_p 2950: the enthalpy is the state,
the temperature read back from the table).

A drop below T-burn exchanges as an inert one and heats; from T-burn it burns on the d^2-law of a drop whose vapour
burns in the gas ([GG09] eqs. 4-5, tests/REFERENCES.md) and every heat it receives feeds the burning, so its
temperature stays where it reached T-burn. With no slip (Re = 0, Ranz-Marshall Nu = 2) the rate is constant:

    d^2(t) = d_a^2 - K (t - t_a),   K = 4 k_g Nu ln(1 + B)/(rho_p c_pg),   B = (c_pg (T_g - T_b) + q-comb)/Lv

(at T_b = T-burn: B = 2.898087, K = 9.549606e-08 m^2/s). A burning drop is removed at the first accepted step whose end
has its mass at or below that of a drop of d-burnout, and what it still holds is handed to the gas.

In each run:
  D1  the solver log holds no failure: no `Run_ODESolver err`, `non-finite state`, `maxIter` or `EXIT OF` line
  D2  25 exit records, IDs 1-25, each removed inside the box: bc 0 and X < Lx = 0.15 m
  D3  each exit diameter d_e at the removal size: D_LO d-burnout <= d_e <= d-burnout (+ its half-ULP)
  D4  each drop reaches T-burn and holds it: the injection row prints T0, every later row and the exit record print one
      and the same T_b, with T-burn - hu(T) <= T_b <= T-burn + TOL_T
  D5  the d^2-law: every burning row after the first and the exit record against the closed form anchored at the
      drop's first burning row (t_a, d_a), K at the drop's own T_b: |d^2 - (d_a^2 - K (t - t_a))| <= hu(d) (2 d +
      hu(d)) + hu(d_a) (2 d_a + hu(d_a)) + K (hu(t) + hu(t_a)) + K (t - t_a) e_T + E_INT (d^2 + d_a^2), hu the
      half-ULP of each printed value (d E13.6, t ES 8 decimals, T F12.6), e_T = c_pg hu(T)/(Lv (1 + B) ln(1 + B)) the
      relative change of K over T's print, E_INT = 1e-9 the integrator's (rtol 1e-11, accumulated over the steps)
  D6  the mass books: the sum of wdot(A) over the 1500 cells of source.tec = the mass flow injected, every stream
      burning out inside the box and its remnant handed: the sum over the 25 inlet cells of col1 rho_g u_g A, col1
      the first column of each 401 record of bc.txt, A = 1e-4 m^2 the inlet face of a cell (G = 0), within TOL_M
      relative (round-off of the source sums: 1500 cells, the mollifier's passes)

Tolerances. D5 is the print floor (the integrator's rtol 1e-11 is below it): measured, the largest residual is 0.966
of it in this run and 0.935 in ref-varcp/, the rounding of a drop's anchor shared by all its rows. TOL_T and D_LO are
not set by a closed form: how far past T-burn the step that reaches it ends, and how far below the removal mass the
step that crosses it ends, depend on the step sizes. Measured: the largest T_b - T-burn 1.96e-2 K, the smallest
d_e/d-burnout 0.9927. TOL_T = 1 K (half the injection's 2 K below T-burn) and D_LO = 0.5 bound them for any step
control, while a drop that is not removed reaches d_e < 1e-11 m or NaN.
"""
import math
import os
import re
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "tools"))
from datcols import parcel_zones, zones  # noqa: E402

IGLOO = os.environ.get("IGLOO_EXE", os.path.abspath(os.path.join(HERE, "..", "..", "..", "bin", "IGLOO")))
VARCP = os.path.join(HERE, "..", "..", "standard", "temp-relax-varcp", "INPUT", "properties.dat")

# ---- inputs (SI) ----
T_G, K_G, GAM, R_G = 600.0, 0.026, 1.4, 287.0
RHO_P = 2950.0
T0, T_BURN, Q_COMB, LV, D_BURN = 300.0, 302.0, 3.1e7, 1.08e7, 3.0e-6
RHO_G, U_G, A_CELL = 1.2, 10.0, 1.0e-4
NU = 2.0
LX = 0.15
N_DROPS = 25

CPG = GAM * R_G / (GAM - 1.0)
B = (CPG * (T_G - T_BURN) + Q_COMB) / LV
K = 4.0 * K_G * NU * math.log(1.0 + B) / (RHO_P * CPG)

TOL_T = 1.0
D_LO = 0.5
TOL_M = 1.0e-12
E_INT = 1.0e-9
FAIL_LINES = ("Run_ODESolver err", "non-finite state", "maxIter", "EXIT OF")


def hu_e(v):
    """Half-ULP of an E13.6 value (0.dddddd E+ee); a printed zero is exact."""
    return 0.0 if v == 0.0 else 0.5e-6 * 10.0 ** (math.floor(math.log10(abs(v))) + 1)


def hu_f6(v):
    """Half-ULP of an F12.6 value."""
    return 0.5e-6


def hu_es8(v):
    """Half-ULP of an ES value with 8 decimals (d.dddddddd E+ee)."""
    return 0.0 if v == 0.0 else 0.5e-8 * 10.0 ** math.floor(math.log10(abs(v)))


def finite(*vals):
    return all(math.isfinite(v) for v in vals)


def source_sum(path):
    """The sum of wdot(A) over the cells of source.tec."""
    lines = open(path).read().splitlines()
    zi = next(i for i, l in enumerate(lines) if l.strip().upper().startswith("ZONE"))
    names = lines[zi - 1].split("=", 1)[1].split()
    I, J, Kk = (int(re.search(r"%s=\s*(\d+)" % a, lines[zi]).group(1)) for a in "IJK")
    nn, nc = I * J * Kk, (I - 1) * (J - 1) * (Kk - 1)
    vals = [float(v) for l in lines[zi + 1:] for v in l.split()]
    k = names.index("wdot(A)") - 3
    w = vals[3 * nn + k * nc: 3 * nn + (k + 1) * nc]
    return sum(w), len(w)


def injected_flow(path):
    """The particle mass flow of the 401 records of bc.txt: col1 (mdot_gas + G A), mdot_gas = rho_g u_g A, G = 0."""
    lines = open(path).read().splitlines()
    cols = [float(lines[n + 1].split()[0]) for n, l in enumerate(lines) if l.split()[-1:] == ["401"]]
    return sum(c * RHO_G * U_G * A_CELL for c in cols), len(cols)


def check_run(d, name, checks):
    """D1-D6 on the run in directory d."""
    log = open(os.path.join(d, "run_out.txt"), errors="replace").read().splitlines()
    bad = [l.strip() for l in log if any(s in l for s in FAIL_LINES)]
    checks.append(("%s D1 solver failures in run_out.txt: %d%s"
                   % (name, len(bad), (" (first: %s)" % bad[0]) if bad else ""), not bad))

    ex = {r["ID"]: r for z in zones(os.path.join(d, "OUTPUT", "outloc-A.dat")) for r in z["rows"]}
    ids_ok = sorted(ex) == list(range(1, N_DROPS + 1))
    inside = [i for i, r in ex.items() if r["bc"] == 0 and r["X"] < LX]
    checks.append(("%s D2 exit records %d, IDs 1-%d %s, removed inside the box (bc 0, X < %.2f m) %d"
                   % (name, len(ex), N_DROPS, "present" if ids_ok else "NOT all", LX, len(inside)),
                   ids_ok and len(inside) == N_DROPS))

    de = {i: r["d<sub>p</sub>"] for i, r in ex.items()}
    off = [i for i, v in de.items() if not (finite(v) and D_LO * D_BURN <= v <= D_BURN + hu_e(D_BURN))]
    fin = [v for v in de.values() if finite(v)]
    checks.append(("%s D3 exit diameter in [%.2f, 1] x d-burnout: %d of %d outside (d_e/d-burnout %s)"
                   % (name, D_LO, len(off), len(de),
                      ("%.6f .. %.6f" % (min(fin) / D_BURN, max(fin) / D_BURN)) if fin else "none finite"),
                   not off and len(de) == N_DROPS))

    trj = {z["id"]: z["rows"] for z in parcel_zones(os.path.join(d, "OUTPUT", "trajectories-A.dat"))}
    t_bad, d_bad, worst_t, worst_r, nrows = [], [], 0.0, 0.0, 0
    for i in sorted(ex):
        rows = trj.get(i, [])
        e = ex[i]
        if len(rows) < 2 or rows[0]["T"] != T0:
            t_bad.append(i)
            continue
        burn = [(r["time"], r["d<sub>p"], r["T"]) for r in rows[1:]] + [(e["time"], e["d<sub>p</sub>"], e["T"])]
        Tb = burn[0][2]
        if not (all(b[2] == Tb for b in burn) and T_BURN - hu_f6(T_BURN) <= Tb <= T_BURN + TOL_T):
            t_bad.append(i)
        worst_t = max(worst_t, Tb - T_BURN) if finite(Tb) else float("nan")
        bb = (CPG * (T_G - Tb) + Q_COMB) / LV
        kb = 4.0 * K_G * NU * math.log(1.0 + bb) / (RHO_P * CPG)
        e_t = CPG * hu_f6(Tb) / (LV * (1.0 + bb) * math.log(1.0 + bb))
        ta, da, _ = burn[0]
        for t, dd, _ in burn[1:]:
            nrows += 1
            r = float("inf")
            if finite(dd, t, da, ta, kb):
                pred = da * da - kb * (t - ta)
                tol = (hu_e(dd) * (2 * dd + hu_e(dd)) + hu_e(da) * (2 * da + hu_e(da)) + kb * (hu_es8(t) + hu_es8(ta))
                       + kb * (t - ta) * e_t + E_INT * (dd * dd + da * da))
                r = abs(dd * dd - pred) / tol
            worst_r = max(worst_r, r)
            if not r <= 1.0:
                d_bad.append(i)
                break
    checks.append(("%s D4 T0 = %.1f K at injection, one T from the first burning row to the exit, T-burn <= T <= "
                   "T-burn + %.2f K: %d of %d drops off (largest T - T-burn %.3e K)"
                   % (name, T0, TOL_T, len(t_bad), len(ex), worst_t), not t_bad))
    checks.append(("%s D5 d^2-law from the first burning row, K at the drop's T_b, %d rows and exit records: %d drops "
                   "off, largest residual/tolerance %.3f" % (name, nrows, len(d_bad), worst_r),
                   not d_bad and nrows > 0))

    w, nc = source_sum(os.path.join(d, "OUTPUT", "source.tec"))
    tot, ninj = injected_flow(os.path.join(d, "INPUT", "bc.txt"))
    rel = abs(w - tot) / tot
    checks.append(("%s D6 sum of wdot(A) over %d cells %.15e kg/s vs the flow of the %d inlet cells %.15e kg/s: "
                   "rel %.2e <= %.0e" % (name, nc, w, ninj, tot, rel, TOL_M), rel <= TOL_M and ninj == N_DROPS))


def sibling(name, table):
    """This case's input run in HERE/<name>/ with another properties table."""
    d = os.path.join(HERE, name)
    shutil.rmtree(d, ignore_errors=True)
    os.makedirs(os.path.join(d, "OUTPUT"))
    os.makedirs(os.path.join(d, "INPUT"))
    for f in ("solfile.tec", "bc.txt", "phase.txt", "models.txt"):
        os.symlink(os.path.realpath(os.path.join(HERE, "INPUT", f)), os.path.join(d, "INPUT", f))
    os.symlink(os.path.realpath(table), os.path.join(d, "INPUT", "properties.dat"))
    shutil.copy(os.path.join(HERE, "input.ini"), os.path.join(d, "input.ini"))
    with open(os.path.join(d, "run_out.txt"), "w") as out, open(os.path.join(d, "run_err.txt"), "w") as err:
        rc = subprocess.run([IGLOO], cwd=d, stdout=out, stderr=err).returncode
    return d, rc


def main():
    checks = []
    print("burnout-direct: B = %.6f, K = %.6e m^2/s, c_pg = %.1f J/kg/K" % (B, K, CPG))
    check_run(HERE, "T-state", checks)
    d, rc = sibling("ref-varcp", VARCP)
    if rc != 0:
        checks.append(("h-state sibling run (ref-varcp) exit %d" % rc, False))
    else:
        log = open(os.path.join(d, "run_out.txt"), errors="replace").read()
        checks.append(("h-state ref-varcp solves the enthalpy (the run's banner)", "Solving enthalpy equation" in log))
        check_run(d, "h-state", checks)

    ok = True
    for text, passed in checks:
        print("[%s] %s" % ("PASS" if passed else "FAIL", text))
        ok = ok and passed
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
