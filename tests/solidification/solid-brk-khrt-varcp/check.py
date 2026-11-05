#!/usr/bin/env python3
"""Independent oracle for Reitz KH-RT breakup on a SOLIDIFYING material (ODE model 7) with a varying liquid heat
capacity and a varying density.

solid-brk-khrt's case with its own INPUT/properties.dat: Cp = 1250 + 0.5 (T - 2400), the Enthalpy column its exact
integral, Density = 2950 + 0.5 (2400 - T). The drops slip (Stokes drag, Ranz-Marshall), so there is no closed form
in x; the gates are solid-brk-khrt's books with the table's enthalpy and density, read from the table file itself.
The state stays the temperature: the KH shed, a child's first state, the density at a breakup event and the
breakup's properties read it as a temperature, and so do the trajectory, the scatter cloud and the Eulerian field.

Gates
-----
K1  solid-brk-khrt's: family 2 sheds no child, family 1 does, the logged children are family 1's child exit records,
    every child with two trajectory rows at least
K2  from a parcel's first frozen row on (parents and children) every row carries the m of that row (the density
    varies, so d follows T at a constant mass); every family-2 parcel's m has decreased on that row (stripping in
    the liquid, before the nucleation); on every trajectory row |m - rho(T) pi d^3/6| within the print bound
K3  solid-brk-khrt's: every exit record at the outlet, each family's exit flows sum to N_INJ injected flows
K4  sum of source.tec E = sum over the injected flows of mdot (h(T_0) + v_0^2/2) minus the sum over the exit
    records of mdot (h_exit + |u_exit|^2/2), h from the table by the phase of the record
K5  solid-brk-khrt's Eulerian mass per family
K6  every scatter row: 0 < m <= the larger injected droplet mass, |m - rho(T) pi d^3/6| within the print bound,
    0 <= f <= 1; rows with 0 < f < 1 and with f = 1 both occur
K7  every Eulerian cell holding particles has T_p between the gas temperature and the hotter injection temperature,
    within 1e-9 relative (the rounding of the cell's weighted mean)
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


khrt = load("solid_brk_khrt_check", os.path.join(HERE, "..", "solid-brk-khrt", "check.py"))
tab = khrt.tab
datcols = tab.datcols
TABLE = os.path.join(HERE, "INPUT", "properties.dat")
HALF_ULP_E13 = khrt.HALF_ULP_E13
TOL_K7 = 1.0e-9           # the rounding of the cell's rho-weighted mean temperature


def read_table(path):
    """{T node: (h, rho)} from the Temperature, Enthalpy and Density columns of the table file."""
    names, out = None, {}
    for line in open(path):
        s = line.strip()
        if s.upper().startswith("VARIABLES"):
            names = re.findall(r'"([^"]*)"', s)
            continue
        if not s or names is None or not (s[0].isdigit() or s[0] in "+-."):
            continue
        v = [float(t) for t in s.split()]
        out[int(round(v[names.index("Temperature")]))] = (v[names.index("Enthalpy")], v[names.index("Density")])
    return out


TAB = read_table(TABLE)
T_LO, T_HI = min(TAB), max(TAB)


def interp(T, k):
    """Linear between the nodes, the end segments extended for h (k = 0), the end values for rho (k = 1)."""
    if k == 1:
        T = min(max(T, float(T_LO)), float(T_HI))
        if T >= T_HI:
            return TAB[T_HI][1]
    i = min(max(int(T), T_LO), T_HI - 1)
    return TAB[i][k] + (TAB[i + 1][k] - TAB[i][k]) * (T - i)


def h_liq(T):
    return interp(T, 0)


def rho(T):
    return interp(T, 1)


def h_mass(T, f):
    """Enthalpy per unit mass at the phase (f, T) gives, on the table."""
    ph = datcols.solid_phase(f, T, tab.T_M)
    if ph == 3:
        return h_liq(tab.T_M) - tab.H_FUS - tab.C_S * (tab.T_M - T)
    if ph == 2:
        return h_liq(tab.T_M) - f * tab.H_FUS
    return h_liq(T)


M_0 = {c: rho(tab.T_0[c]) * math.pi / 6.0 * tab.D_0 ** 3 for c in tab.T_0}


def md_ratio(d_tok, m_tok, T):
    """|m - rho(T) pi d^3/6| over its print bound (E13.6 half-ULPs of d and m)."""
    d, m = float(d_tok), float(m_tok)
    return abs(m - rho(T) * math.pi / 6.0 * d ** 3) / m / (3.0 * khrt.ulp_rel(d_tok) + khrt.ulp_rel(m_tok))


def check_mass_rows(fails):
    """K2's m-d relation on every trajectory row, from the printed tokens: (rows, worst ratio)."""
    n, worst, nbad = 0, 0.0, 0
    with open(tab.TRAJ) as fh:
        names = None
        for line in fh:
            s = line.strip()
            low = s.lower()
            if low.startswith("variables"):
                names = re.findall(r'"([^"]*)"', s)
                iD, iM, iT = names.index("d<sub>p"), names.index("m<sub>p"), names.index("T")
                continue
            if not s or low.startswith("zone"):
                continue
            c = s.split()
            n += 1
            r = md_ratio(c[iD], c[iM], float(c[iT]))
            worst = max(worst, r)
            if r > 1.0:
                nbad += 1
                if nbad <= 3:
                    fails.append(f"K2 trajectory row {n}: d={c[iD]} m={c[iM]} T={c[iT]}: m - rho(T) pi d^3/6 at {r:.3g} "
                                 f"of the print bound")
    if nbad:
        fails.append(f"K2 {nbad} of {n} trajectory rows violate the m-d relation")
    return n, worst


def check_frozen_m(traj, fails):
    """K2: from the first frozen row on every row carries that row's m; family 2's m decreased on that row."""
    n_rows, unstripped = 0, []
    for (fam, pid), rows in sorted(traj.items()):
        k = tab.first_frozen(rows)
        if k is None:
            fails.append(f"K2 family {fam} ID={pid}: no frozen row (phase >= 2) on its {len(rows)} rows")
            continue
        mk = rows[k]["m<sub>p"]
        bad = [r for r in rows[k + 1:] if r["m<sub>p"] != mk]
        n_rows += len(rows) - k
        if bad:
            fails.append(f"K2 family {fam} ID={pid}: {len(bad)} rows after the first frozen row (t={rows[k]['t']:.6e}, "
                         f"m={mk:.6e}) change m, first t={bad[0]['t']:.6e} m={bad[0]['m<sub>p']:.6e}")
        if fam == 2 and (k == 0 or not rows[k]["m<sub>p"] < rows[k - 1]["m<sub>p"]):
            unstripped.append(pid)
    if unstripped:
        fails.append(f"K2 non-vacuity: {len(unstripped)} family-2 parcels with no m decrease on the first frozen row, "
                     f"first ID={unstripped[0]}")
    return n_rows, len(unstripped)


def check_energy(exits, fails):
    """K4: the global energy telescoping on the table's enthalpy."""
    wdot, _, _, _, E = tab.read_source(tab.SRC)
    e_in = sum(tab.N_INJ * tab.MDOT[c] * (h_liq(tab.T_0[c]) + 0.5 * tab.V_0 * tab.V_0) for c in tab.MDOT)
    e_out = sum(r["mdot"] * (h_mass(r["T"], r["f"]) + 0.5 * r["|u<sub>p</sub>|"] ** 2)
                for rr in exits.values() for r in rr)
    ref, got = e_in - e_out, sum(E)
    rel = abs(got - ref) / abs(ref) if ref else float("inf")
    if rel > tab.TOL_E:
        fails.append(f"K4 sum(E) = {got:.9e} W vs the closed form {ref:.9e} W (rel {rel:.3e} > {tab.TOL_E:.0e})")
    if any(w != 0.0 for w in wdot):
        fails.append("K4 wdot /= 0: a constant-mass material exchanges no mass")
    return got, ref, rel


def check_scatter(fails):
    """K6, streamed: (rows, worst ratio to the print bound, largest m, rows with 0 < f < 1, rows with f = 1)."""
    n, worst, mmax, nbad, nmid, nsol = 0, 0.0, 0.0, 0, 0, 0
    m_top = max(M_0.values()) * (1.0 + 2.0 * HALF_ULP_E13)
    with open(khrt.SCAT) as fh:
        names = None
        for line in fh:
            s = line.strip()
            low = s.lower()
            if low.startswith("variables"):
                names = re.findall(r'"([^"]*)"', s)
                iD, iM, iF, iT = names.index("d<sub>p"), names.index("m<sub>p"), names.index("f"), names.index("T")
                continue
            if not s or low.startswith("zone"):
                continue
            c = s.split()
            if names is None or len(c) != len(names):
                raise ValueError(f"{khrt.SCAT}: a row of {len(c)} values for {len(names or [])} names")
            d, m, f = float(c[iD]), float(c[iM]), float(c[iF])
            n += 1
            nmid += 0.0 < f < 1.0
            nsol += f == 1.0
            ratio = float("inf")
            if d > 0.0 and 0.0 < m <= m_top and 0.0 <= f <= 1.0:
                ratio = md_ratio(c[iD], c[iM], float(c[iT]))
                mmax = max(mmax, m)
            worst = max(worst, ratio)
            if ratio > 1.0:
                nbad += 1
                if nbad <= 3:
                    fails.append(f"K6 scatter row {n}: d={c[iD]} m={c[iM]} T={c[iT]} f={c[iF]} (m - rho(T) pi d^3/6 at "
                                 f"{ratio:.3g} of the print bound, or m, f out of range)")
    if nbad:
        fails.append(f"K6 {nbad} of {n} scatter rows violate the bounds")
    if n == 0:
        fails.append("K6 non-vacuity: no scatter row")
    if nmid == 0 or nsol == 0:
        fails.append(f"K6 non-vacuity: {nmid} rows with 0 < f < 1, {nsol} with f = 1")
    return n, worst, mmax, nmid, nsol


def euler_temperature(path):
    """(cells holding particles, lowest T_p, highest T_p) from the cell-centred block of an euler<c>.tec."""
    lines = open(path).read().splitlines()
    zone = next(l for l in lines if l.strip().startswith("ZONE"))
    ni, nj, nk = (int(re.search(r"%s=(\d+)" % c, zone).group(1)) for c in "IJK")
    nnod, ncel = ni * nj * nk, (ni - 1) * (nj - 1) * (nk - 1)
    vals = [float(v) for v in " ".join(lines[lines.index(zone) + 1:]).split()]
    rho_p = vals[3 * nnod:3 * nnod + ncel]
    t_p = vals[3 * nnod + 4 * ncel:3 * nnod + 5 * ncel]
    held = [t for r, t in zip(rho_p, t_p) if r > 0.0]
    return len(held), min(held, default=float("nan")), max(held, default=float("nan"))


def main():
    try:
        traj, exits = tab.parcels(tab.TRAJ), tab.parcels(tab.EXIT)
        log = open(tab.LOG, errors="replace").read()
    except (FileNotFoundError, StopIteration, ValueError, KeyError) as e:
        print(f"[FAIL] cannot read the outputs: {e}")
        return 1
    fails = []

    logged = khrt.children_logged(log)
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
    few_rows = [pid for pid in kids[1] if len(traj.get((1, pid), [])) < 2]
    if few_rows:
        fails.append(f"K1 {len(few_rows)} of {len(kids[1])} children with fewer than two trajectory rows, "
                     f"first ID={few_rows[0]}")

    n2, n_unstripped = check_frozen_m(traj, fails)
    try:
        n2r, w2r = check_mass_rows(fails)
    except (FileNotFoundError, ValueError) as e:
        fails.append(f"K2 {e}")
        n2r, w2r = 0, float("inf")
    w3 = tab.check_flows(exits, fails, "K3", per_record=False)
    got, ref, rel4 = check_energy(exits, fails)

    k5, k7 = {}, {}
    for c in tab.MDOT:
        try:
            em = khrt.euler_mass(khrt.EUL.format(c))
            k7[c] = euler_temperature(khrt.EUL.format(c))
        except (FileNotFoundError, StopIteration, ValueError) as e:
            fails.append(f"K5 family {c}: {e}")
            continue
        ft = sum(r["mdot"] * r["t"] for (f, _), rr in exits.items() if f == c for r in rr)
        k5[c] = (em, ft, abs(em - ft) / ft if ft > 0.0 else float("inf"))
        if k5[c][2] > tab.TOL_E:
            fails.append(f"K5 family {c}: sum rho_p V = {em:.9e} kg vs sum mdot t = {ft:.9e} kg (rel {k5[c][2]:.3e})")
        n7, lo, hi = k7[c]
        if n7 == 0 or not (tab.T_G * (1.0 - TOL_K7) <= lo and hi <= max(tab.T_0.values()) * (1.0 + TOL_K7)):
            fails.append(f"K7 family {c}: {n7} cells, T_p in [{lo:.3f}, {hi:.3f}] K, not within [{tab.T_G}, "
                         f"{max(tab.T_0.values())}]")

    try:
        n6, w6, mmax, nmid, nsol = check_scatter(fails)
    except (FileNotFoundError, ValueError) as e:
        fails.append(f"K6 {e}")
        n6, w6, mmax, nmid, nsol = 0, float("inf"), 0.0, 0, 0

    good = tab.check_injection(traj, exits, fails, "K-guard")
    for w in ("solidification --> supercooling + recalescence (model 7)", "breakup     --> Reitz-KHRT",
              "eulerian field included"):
        if log.count(w) != 2:
            fails.append(f"K-guard witness: the log does not report '{w}' for both groups")
    if "Solving enthalpy equation" in log:
        fails.append("K-guard witness: the log reports an enthalpy state for a solidifying material")

    print(f"oracle: injected droplet mass {M_0[1]:.6e} / {M_0[2]:.6e} kg (rho {rho(tab.T_0[1]):.3f} / "
          f"{rho(tab.T_0[2]):.3f}), h(T_0) {h_liq(tab.T_0[1]):.3f} / {h_liq(tab.T_0[2]):.3f} J/kg")
    print(f"K1 children logged per family {dict((g, sum(v)) for g, v in logged.items())}, child exit records "
          f"{len(kids[1])} / {len(kids[2])}")
    print(f"K2 {n2} rows from the first frozen row on keep their m; {25 - n_unstripped} of 25 family-2 parcels "
          f"stripped before the nucleation; {n2r} trajectory rows, m-d worst at {w2r:.3f} of the print bound")
    print(f"K3 {sum(len(v) for v in exits.values())} exit records; worst family flow deviation {w3:.3e} "
          f"(tol {tab.TOL_E:.0e})")
    print(f"K4 sum(E) = {got:.9e} W, closed form {ref:.9e} W, rel {rel4:.3e}")
    for c, (em, ft, rel) in sorted(k5.items()):
        print(f"K5 family {c}: sum rho_p V = {em:.9e} kg, sum mdot t = {ft:.9e} kg, rel {rel:.3e}")
    for c, (n7, lo, hi) in sorted(k7.items()):
        print(f"K7 family {c}: {n7} cells, T_p in [{lo:.3f}, {hi:.3f}] K")
    print(f"K6 {n6} scatter rows, worst at {w6:.3f} of the print bound, largest m {mmax:.6e} kg; {nmid} rows with "
          f"0 < f < 1, {nsol} with f = 1")
    print(f"guard: parcels verified per family {good}")
    for f in fails[:40]:
        print(f"[FAIL] {f}")
    if fails:
        print(f"\n[FAIL] {len(fails)} violation(s)")
        return 1
    print("\n[PASS] KH-RT on the liquid only with the table's heat capacity and density: the frozen drops keep their "
          "mass; flows, energy, the Eulerian mass and temperature balance.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
