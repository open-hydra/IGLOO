#!/usr/bin/env python3
"""Independent oracle for Reitz KH-RT breakup on a SOLIDIFYING material (ODE model 7, with KH children).

solid-brk-tab's box, inlet and material (its check.py supplies the inputs and the shared gates) under KH-RT:
family 1 is molten, sheds KH children and breaks up by RT before it nucleates; family 2 nucleates within a
few microseconds of the injection and must shed nothing. KH stripping acts in the liquid only (phases 0 and
1), so it changes a parcel's d up to its first frozen row k (the first row with phase >= 2) and never after;
a child is born in its parent's phase -- here always the liquid above T_m -- and integrated from there (its
first written row lies at the next face, already solid).

Gates
-----
K1  family 2 sheds no child (its group reports no "number of children"), family 1 does, and the children the
    log reports are the child IDs (> N_INJ) of family 1's exit records, one record each, every child with
    trajectory rows
K2  every row from a parcel's first frozen row on carries the d of that row, parents and children; every
    family-2 parcel's d has changed on that row (stripping in the liquid, before the nucleation)
K3  every exit record ends at the outlet (code 400), and each family's exit flows (parents and children) sum
    to N_INJ injected flows
K4  the energy balance of solid-brk-tab's T5 over parents and children
K5  the Eulerian mass of each family: sum over cells of rho_p V (euler<c>.tec) = sum over its exit records of
    mdot t (mollify off; the flow changes only at a KH shed, where the child takes what the parent gives)
K6  every scatter row: 0 < d <= d_0, |m - rho pi d^3/6| within the print bound of the row (E13.6 half-ULPs of d
    and m), 0 <= f <= 1; rows with 0 < f < 1 and rows with f = 1 both occur (model 7 writes its f). The bound
    d <= d_0 rests on the constant density (no growth path here); the m-d relation witnesses the branch that
    derives d from m, it is not a physical oracle
"""
import importlib.util
import math
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


tab = load("solid_brk_tab_check", os.path.join(HERE, "..", "solid-brk-tab", "check.py"))
datcols = tab.datcols
SCAT, EUL = "OUTPUT/scatter-A.dat", "OUTPUT/euler{}.tec"
HALF_ULP_E13 = 0.5e-6     # E13.6: 0.dddddd mantissa


def euler_mass(path):
    """Sum over cells of rho_p times the cell volume, from the cell-centred block of an euler<c>.tec."""
    lines = open(path).read().splitlines()
    zone = next(l for l in lines if l.strip().startswith("ZONE"))
    ni, nj, nk = (int(re.search(r"%s=(\d+)" % c, zone).group(1)) for c in "IJK")
    nnod, ncel = ni * nj * nk, (ni - 1) * (nj - 1) * (nk - 1)
    vals = [float(v) for v in " ".join(lines[lines.index(zone) + 1:]).split()]
    x, y, z = vals[:nnod], vals[nnod:2 * nnod], vals[2 * nnod:3 * nnod]
    rho = vals[3 * nnod:3 * nnod + ncel]
    tot = 0.0
    for k in range(nk - 1):
        for j in range(nj - 1):
            for i in range(ni - 1):
                n0 = i + ni * (j + nj * k)
                vol = (x[n0 + 1] - x[n0]) * (y[n0 + ni] - y[n0]) * (z[n0 + ni * nj] - z[n0])
                tot += rho[i + (ni - 1) * (j + (nj - 1) * k)] * vol
    return tot


def children_logged(log):
    """{group: [children counts of its loops]} from the group blocks of the log."""
    out, g = {}, None
    for line in log.splitlines():
        m = re.search(r"Group\s+(\d+)\s+=> number of particles", line)
        if m:
            g = int(m.group(1))
            out[g] = []
            continue
        m = re.search(r"number of children =\s*(\d+)", line)
        if m and g is not None:
            out[g].append(int(m.group(1)))
    return out


def ulp_rel(tok):
    """The E13.6 half-ULP of a printed token, relative to its value."""
    return HALF_ULP_E13 * 10.0 ** int(tok[-3:]) / abs(float(tok))


def check_scatter(fails):
    """K6, streamed: (rows, worst ratio to the print bound, largest d, rows with 0 < f < 1, rows with f = 1)."""
    n, worst, dmax, nbad, nmid, nsol = 0, 0.0, 0.0, 0, 0, 0
    with open(SCAT) as fh:
        names = None
        for line in fh:
            s = line.strip()
            low = s.lower()
            if low.startswith("variables"):
                names = re.findall(r'"([^"]*)"', s)
                iD, iM, iF = names.index("d<sub>p"), names.index("m<sub>p"), names.index("f")
                continue
            if not s or low.startswith("zone"):
                continue
            c = s.split()
            if names is None or len(c) != len(names):
                raise ValueError(f"{SCAT}: a row of {len(c)} values for {len(names or [])} names")
            d, m, f = float(c[iD]), float(c[iM]), float(c[iF])
            n += 1
            nmid += 0.0 < f < 1.0
            nsol += f == 1.0
            ok = 0.0 < d <= tab.D_0 and m > 0.0 and 0.0 <= f <= 1.0
            ratio = float("inf")
            if ok:
                ratio = abs(m - tab.RHO_P * math.pi / 6.0 * d ** 3) / m / (3.0 * ulp_rel(c[iD]) + ulp_rel(c[iM]))
                dmax = max(dmax, d)
            worst = max(worst, ratio)
            if ratio > 1.0:
                nbad += 1
                if nbad <= 3:
                    fails.append(f"K6 scatter row {n}: d={c[iD]} m={c[iM]} f={c[iF]} (m - rho pi d^3/6 at {ratio:.3g} of "
                                 f"the print bound, or d, m, f out of range)")
    if nbad:
        fails.append(f"K6 {nbad} of {n} scatter rows violate the bounds")
    if n == 0:
        fails.append("K6 non-vacuity: no scatter row")
    if nmid == 0 or nsol == 0:
        fails.append(f"K6 non-vacuity: {nmid} rows with 0 < f < 1, {nsol} with f = 1")
    return n, worst, dmax, nmid, nsol


def main():
    try:
        traj, exits = tab.parcels(tab.TRAJ), tab.parcels(tab.EXIT)
        log = open(tab.LOG, errors="replace").read()
    except (FileNotFoundError, StopIteration, ValueError, KeyError) as e:
        print(f"[FAIL] cannot read the outputs: {e}")
        return 1
    fails = []

    # K1: children by family, from the log and the exit records
    logged = children_logged(log)
    kids = {c: sorted(pid for (f, pid) in exits if f == c and pid > tab.N_INJ) for c in tab.MDOT}
    multi = [k for k, v in exits.items() if len(v) != 1]
    if sorted(logged) != [1, 2]:
        fails.append(f"K1 the log has group blocks {sorted(logged)}, not [1, 2]")
    if logged.get(2):
        fails.append(f"K1 family 2 sheds children: {logged[2]}")
    if not logged.get(1) or sum(logged[1]) <= 0:
        fails.append("K1 family 1 sheds no child")
    if sum(logged.get(1, [])) != len(kids[1]) or kids[2]:
        fails.append(f"K1 children logged {sum(logged.get(1, []))}, child exit records {len(kids[1])} (family 1) and "
                     f"{len(kids[2])} (family 2)")
    if multi:
        fails.append(f"K1 {len(multi)} parcels with other than one exit record, first {multi[0]}")
    no_rows = [pid for pid in kids[1] if (1, pid) not in traj]
    if no_rows:
        fails.append(f"K1 {len(no_rows)} of {len(kids[1])} children have no trajectory row, first ID={no_rows[0]}")

    # K2: the frozen d, parents and children; family 2 stripped before its nucleation
    n2, on_k, _ = tab.check_frozen_d(traj, fails, "K2")
    unstripped, n_fam2 = [], 0
    for (fam, pid), rows in traj.items():
        if fam == 2:
            n_fam2 += 1
            k = tab.first_frozen(rows)
            if k is None or k == 0 or rows[k]["d<sub>p"] >= rows[k - 1]["d<sub>p"]:
                unstripped.append(pid)
    if unstripped:
        fails.append(f"K2 non-vacuity: {len(unstripped)} family-2 parcels with no d decrease on the first frozen row, "
                     f"first ID={unstripped[0]}")

    w3 = tab.check_flows(exits, fails, "K3", per_record=False)
    got, ref, rel4 = tab.check_energy(exits, fails, "K4")

    # K5: the Eulerian mass per family
    k5 = {}
    for c in tab.MDOT:
        try:
            em = euler_mass(EUL.format(c))
        except (FileNotFoundError, StopIteration, ValueError) as e:
            fails.append(f"K5 family {c}: {e}")
            continue
        ft = sum(r["mdot"] * r["t"] for (f, _), rr in exits.items() if f == c for r in rr)
        k5[c] = (em, ft, abs(em - ft) / ft if ft > 0.0 else float("inf"))
        if k5[c][2] > tab.TOL_E:
            fails.append(f"K5 family {c}: sum rho_p V = {em:.9e} kg vs sum mdot t = {ft:.9e} kg (rel {k5[c][2]:.3e})")

    try:
        n6, w6, dmax, nmid, nsol = check_scatter(fails)
    except (FileNotFoundError, ValueError) as e:
        fails.append(f"K6 {e}")
        n6, w6, dmax, nmid, nsol = 0, float("inf"), 0.0, 0, 0

    good = tab.check_injection(traj, exits, fails, "K-guard")
    for w in ("solidification --> supercooling + recalescence (model 7)", "breakup     --> Reitz-KHRT",
              "eulerian field included"):
        if log.count(w) != 2:
            fails.append(f"K-guard witness: the log does not report '{w}' for both groups")

    print(f"K1 children logged per family {dict((g, sum(v)) for g, v in logged.items())}, child exit records "
          f"{len(kids[1])} / {len(kids[2])}")
    print(f"K2 {n2} rows from the first frozen row on (parents and children); {on_k} parcels change d on that row; "
          f"{n_fam2 - len(unstripped)} of {n_fam2} family-2 parcels stripped before the nucleation")
    print(f"K3 {sum(len(v) for v in exits.values())} exit records; worst family flow deviation {w3:.3e} "
          f"(tol {tab.TOL_E:.0e})")
    print(f"K4 sum(E) = {got:.9e} W, closed form {ref:.9e} W, rel {rel4:.3e}")
    for c, (em, ft, rel) in sorted(k5.items()):
        print(f"K5 family {c}: sum rho_p V = {em:.9e} kg, sum mdot t = {ft:.9e} kg, rel {rel:.3e}")
    print(f"K6 {n6} scatter rows, worst at {w6:.3f} of the print bound, largest d {dmax:.6e} m; {nmid} rows with "
          f"0 < f < 1, {nsol} with f = 1")
    print(f"guard: parcels verified per family {good}")
    for f in fails[:40]:
        print(f"[FAIL] {f}")
    if fails:
        print(f"\n[FAIL] {len(fails)} violation(s)")
        return 1
    print("\n[PASS] KH-RT acts on the liquid only: the molten family sheds and breaks up, the frozen drops keep their "
          "size; flows, energy and the Eulerian mass balance.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
