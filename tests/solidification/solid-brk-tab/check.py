#!/usr/bin/env python3
"""Independent oracle for TAB breakup on a SOLIDIFYING material (ODE model 6, the breakup event gated by the phase).

Every constant below is a case input (bc.txt, properties.dat, the phase-file tokens, the box gas, input.ini);
nothing physical is read back from production output to build the reference. The particle files are read
by the names of their headers (tools/datcols.py); a parcel's rows are taken in time order (the t column).

Physics
-------
Two inlet families of one material enter the box at 35 m/s against the 10 m/s gas. Family 1 (T_0 = 4 T_g)
breaks up while liquid, then nucleates; family 2 (T_0 = 3.10433 T_g, 1 K above T_nuc) nucleates within a
few microseconds, long before the TAB oscillator reaches breakup, so it never breaks up. Breakup acts on the
liquid only (phases 0 and 1): from the parcel's first frozen row k on -- the first row whose phase,
datcols.solid_phase(f, T, T_m), is 2 (plateau) or 3 (solid) -- its diameter is fixed. A breakup inside the
nucleation cell, before the nucleation, lands on row k itself.
Enthalpy per unit mass (relative table, datum 0): c_l T (liquid, phases 0 and 1), c_l T_m - f h_fus (plateau),
c_l T_m - h_fus - c_s (T_m - T) (solid).

Gates
-----
T1  every family-1 parcel's d decreases at least once, on a row up to its first frozen row k
T2  every family-2 row carries the injected d (the printed value, exactly)
T3  every row from a parcel's first frozen row on carries the d of that row (the printed value, exactly)
T4  every exit record ends at the outlet (code 400) with its family's injected flow, and each family's total
    is N_INJ of them
T5  sum of source.tec E over all cells = sum over the injected flows of mdot (h_0 + v_0^2/2) minus the sum
    over the exit records of mdot (h_exit + |u_exit|^2/2), h by the phase of the exit record; no mass exchange
T6  >= 20 parcels per family with the injection row, a frozen row and one exit record; the model-6 witness
"""
import math
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "tools"))
import datcols                                  # noqa: E402

TRAJ, EXIT, SRC, LOG = "OUTPUT/trajectories-A.dat", "OUTPUT/outloc-A.dat", "OUTPUT/source.tec", "run_out.txt"

# ---- known inputs (SI), NOT read from production output ----------------------------------
T_G, U_G, RHO_G, AREA = 600.0, 10.0, 1.2, 1.0e-4    # box gas (make_box_case.py), inlet cell area
RHO_P, C_L = 2950.0, 1250.0                         # common/properties.dat (Density, Cp)
H_FUS, C_S = 1.07e6, 600.0                          # phase-file tokens h-fus, cp-solid
T_M = 2327.0                                        # T-melt default
D_0 = 30.0e-6                                       # bc.txt rp = 1.5e-5 (both families)
V_0 = 3.5 * U_G                                     # bc.txt kV = 3.5
KT = {1: 4.0, 2: 3.10433}                           # bc.txt kT per family
KRHO = {1: 0.34, 2: 0.17}                           # bc.txt krho per family
N_INJ = 25                                          # inlet cells per family (5 x 5, fsample 1)
BC_OUT = 400                                        # the outlet face's code (bc.txt)
T_0 = {c: KT[c] * T_G for c in KT}
# 401 inlet: mdot = krho / (1 - sum krho) * rho_g u_g A
MDOT = {c: KRHO[c] / (1.0 - sum(KRHO.values())) * RHO_G * U_G * AREA for c in KRHO}

# ---- error model -----------------------------------------------------------------------------
TOL_E = 1.0e-5          # flows and energy (outloc mdot E13.6: half-ULP <= 5e-6 per record; source.tec E23.15)
TOL_T0 = 1.0e-3         # the injection row's T (F12.6)
N_GOOD = 20


def h_mass(T, f):
    """Enthalpy per unit mass at the phase (f, T) gives."""
    ph = datcols.solid_phase(f, T, T_M)
    if ph == 3:
        return C_L * T_M - H_FUS - C_S * (T_M - T)
    if ph == 2:
        return C_L * T_M - f * H_FUS
    return C_L * T


def parcels(path):
    """{(family, ID): [row, ...]} in time order (stable on the file order), family = the zone's group number."""
    out = {}
    for z in datcols.zones(path):
        m = re.search(r"Group\s+(\d+)", z["title"])
        if not m:
            raise ValueError(f"{path}: zone '{z['title']}' names no group")
        fam = int(m.group(1))
        for r in z["rows"]:
            out.setdefault((fam, r["ID"]), []).append(r)
    return {k: sorted(v, key=lambda r: r["t"]) for k, v in out.items()}


def first_frozen(rows):
    """Index of the first row with phase >= 2, None if the parcel never freezes."""
    return next((i for i, r in enumerate(rows) if datcols.solid_phase(r["f"], r["T"], T_M) >= 2), None)


def read_source(path):
    lines = open(path).read().splitlines()
    zi = next(i for i, l in enumerate(lines) if l.strip().upper().startswith("ZONE"))
    I, J, K = (int(re.search(rf"{a}=\s*(\d+)", lines[zi]).group(1)) for a in "IJK")
    vals = [float(t) for l in lines[zi + 1:] for t in l.split()]
    nn, nc = I * J * K, (I - 1) * (J - 1) * (K - 1)
    return [vals[3 * nn + a * nc: 3 * nn + (a + 1) * nc] for a in range(5)]   # wdot Fx Fy Fz E


def check_frozen_d(traj, fails, tag):
    """T3 / K2: from the first frozen row on, every row carries that row's d. Returns (rows checked, parcels whose
    d changed ON the first frozen row, parcels that never froze)."""
    n_rows, on_k, never = 0, 0, 0
    for (fam, pid), rows in sorted(traj.items()):
        k = first_frozen(rows)
        if k is None:
            never += 1
            fails.append(f"{tag} family {fam} ID={pid}: no frozen row (phase >= 2) on its {len(rows)} rows")
            continue
        if k > 0 and rows[k]["d<sub>p"] != rows[k - 1]["d<sub>p"]:
            on_k += 1
        dk = rows[k]["d<sub>p"]
        bad = [r for r in rows[k + 1:] if r["d<sub>p"] != dk]
        n_rows += len(rows) - k
        if bad:
            fails.append(f"{tag} family {fam} ID={pid}: {len(bad)} rows after the first frozen row (t={rows[k]['t']:.6e},"
                         f" d={dk:.6e}) change d, first t={bad[0]['t']:.6e} d={bad[0]['d<sub>p']:.6e}")
    return n_rows, on_k, never


def check_flows(exits, fails, tag, per_record):
    """T4 / K3: every exit record ends at the outlet; each family's exit flows sum to N_INJ injected flows
    (per_record: every record carries one)."""
    worst = 0.0
    off = [(f, pid, r["bc"]) for (f, pid), rr in sorted(exits.items()) for r in rr if r["bc"] != BC_OUT]
    if off:
        fails.append(f"{tag} {len(off)} exit records not at the outlet (code {BC_OUT}), first family {off[0][0]} "
                     f"ID={off[0][1]} code {off[0][2]}")
    for fam in sorted(MDOT):
        recs = [r for (f, _), rr in exits.items() if f == fam for r in rr]
        tot = sum(r["mdot"] for r in recs)
        rel = abs(tot - N_INJ * MDOT[fam]) / (N_INJ * MDOT[fam])
        worst = max(worst, rel)
        if rel > TOL_E:
            fails.append(f"{tag} family {fam}: sum of exit mdot {tot:.9e} vs {N_INJ} x {MDOT[fam]:.9e} (rel {rel:.3e})")
        if per_record:
            for r in recs:
                e = abs(r["mdot"] - MDOT[fam]) / MDOT[fam]
                worst = max(worst, e)
                if e > TOL_E:
                    fails.append(f"{tag} family {fam} ID={r['ID']}: exit mdot {r['mdot']:.6e} vs injected "
                                 f"{MDOT[fam]:.9e} (rel {e:.3e})")
    return worst


def check_energy(exits, fails, tag):
    """T5 / K4: the global energy telescoping over the injected flows and every exit record."""
    wdot, _, _, _, E = read_source(SRC)
    e_in = sum(N_INJ * MDOT[c] * (C_L * T_0[c] + 0.5 * V_0 * V_0) for c in MDOT)
    e_out = sum(r["mdot"] * (h_mass(r["T"], r["f"]) + 0.5 * r["|u<sub>p</sub>|"] ** 2)
                for rr in exits.values() for r in rr)
    ref, got = e_in - e_out, sum(E)
    rel = abs(got - ref) / abs(ref) if ref else float("inf")
    if rel > TOL_E:
        fails.append(f"{tag} sum(E) = {got:.9e} W vs the closed form {ref:.9e} W (rel {rel:.3e} > {TOL_E:.0e})")
    if any(w != 0.0 for w in wdot):
        fails.append(f"{tag} wdot /= 0: a constant-mass material exchanges no mass")
    return got, ref, rel


def check_injection(traj, exits, fails, tag):
    """T6 / non-vacuity: parcels per family with the injection row, a frozen row and one exit record."""
    good = {c: 0 for c in MDOT}
    for (fam, pid), rows in traj.items():
        if pid > N_INJ:
            continue
        r0 = rows[0]
        if r0["t"] != 0.0 or abs(r0["T"] - T_0[fam]) > TOL_T0 or r0["d<sub>p"] != D_0:
            fails.append(f"{tag} family {fam} ID={pid}: first row t={r0['t']}, T={r0['T']}, d={r0['d<sub>p']} is not "
                         f"the injection state (0, {T_0[fam]}, {D_0})")
            continue
        if first_frozen(rows) is not None and len(exits.get((fam, pid), [])) == 1:
            good[fam] += 1
    for c, n in good.items():
        if n < N_GOOD:
            fails.append(f"{tag} family {c}: {n} parcels with the injection row, a frozen row and one exit record "
                         f"(need >= {N_GOOD})")
    return good


def main():
    try:
        traj, exits = parcels(TRAJ), parcels(EXIT)
        log = open(LOG, errors="replace").read()
    except (FileNotFoundError, StopIteration, ValueError, KeyError) as e:
        print(f"[FAIL] cannot read the outputs: {e}")
        return 1
    print(f"oracle: T_0 = {T_0[1]:.3f} / {T_0[2]:.3f} K, injected mdot {MDOT[1]:.9e} / {MDOT[2]:.9e} kg/s, "
          f"d_0 = {D_0:.6e} m, v_0 = {V_0} m/s")
    fails = []

    # T1: family 1 breaks up while liquid
    n1, n1_on_k = 0, 0
    for (fam, pid), rows in sorted(traj.items()):
        if fam != 1:
            continue
        k = first_frozen(rows)
        dec = [i for i in range(1, len(rows)) if rows[i]["d<sub>p"] < rows[i - 1]["d<sub>p"]]
        if k is None or not any(i <= k for i in dec):
            fails.append(f"T1 family 1 ID={pid}: no d decrease up to its first frozen row (k={k}, decreases at {dec[:4]})")
            continue
        n1 += 1
        n1_on_k += k in dec

    # T2: family 2 never breaks up
    bad2 = [(pid, r["t"], r["d<sub>p"]) for (fam, pid), rows in traj.items() if fam == 2 for r in rows
            if r["d<sub>p"] != D_0]
    n2 = sum(len(rows) for (fam, _), rows in traj.items() if fam == 2)
    if bad2:
        fails.append(f"T2 {len(bad2)} family-2 rows with d != {D_0:.6e}, first ID={bad2[0][0]} t={bad2[0][1]:.6e} "
                     f"d={bad2[0][2]:.6e}")

    n3, on_k, _ = check_frozen_d(traj, fails, "T3")
    w4 = check_flows(exits, fails, "T4", per_record=True)
    got, ref, rel5 = check_energy(exits, fails, "T5")
    good = check_injection(traj, exits, fails, "T6")
    if log.count("solidification --> supercooling + recalescence (model 6)") != 2:
        fails.append("T6 witness: the log does not report the model-6 ODE system for both groups")

    print(f"T1 {n1} family-1 parcels break up while liquid ({n1_on_k} of them on the first frozen row itself)")
    print(f"T2 {n2} family-2 rows, {len(bad2)} with d != d_0")
    print(f"T3 {n3} rows from the first frozen row on; {on_k} parcels change d on that row")
    print(f"T4 {sum(len(v) for v in exits.values())} exit records; worst exit-flow deviation {w4:.3e} (tol {TOL_E:.0e})")
    print(f"T5 sum(E) = {got:.9e} W, closed form {ref:.9e} W, rel {rel5:.3e}")
    print(f"T6 parcels verified per family: {good}")
    for f in fails[:40]:
        print(f"[FAIL] {f}")
    if fails:
        print(f"\n[FAIL] {len(fails)} violation(s)")
        return 1
    print("\n[PASS] TAB breaks the molten family up while liquid and never touches a frozen drop; flows and energy "
          "balance.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
