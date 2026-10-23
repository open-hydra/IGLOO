#!/usr/bin/env python3
"""Independent oracle for TAB breakup on a solidifying material with rho-solid (ODE model 6).

solid-brk-tab's box, inlet families and material (its check.py supplies the inputs and the shared gates) with
rho-solid = 3660 kg/m3 on the material line. A drop's mass changes only at a breakup event, and breakup acts on
the liquid only: from a parcel's first frozen row k on (the first row whose phase, datcols.solid_phase(f, T, T_m),
is 2 or 3) its mass is fixed, and on every row its diameter is that of its liquid part at rho_l and its frozen
fraction at rho_s, d = (6 m (a + b f)/pi)^{1/3}, a = 1/rho_l, b = 1/rho_s - 1/rho_l. Family 2 (injected 1 K above
T_nuc) never breaks up: its mass is the injected one on every row, while its diameter shrinks as it freezes.

Gates
-----
W1  every row from a parcel's first frozen row on carries the m of that row (the printed value, exactly); every
    family-2 row carries the injected mass (E13.6 half-ULP)
W2  every row: d = (6 m (a + b f)/pi)^{1/3} from the row's m and f, within the E13.6 half-ULP of d plus a third
    of the E13.6 relative half-ULP of m plus |dln d/df| times the F10.6 half-ULP of f; >= 20 rows per family with
    0 < f < 1 and with f = 1 on which d differs from the liquid-density diameter by more than that budget
W3  every family-1 parcel's m decreases at least once, on a row up to its first frozen row k (breakup while
    liquid; d alone cannot witness it, since freezing shrinks d on row k)
T4, T5, T6  solid-brk-tab's: exit flows, the energy balance, the injection state and the model-6 witnesses
"""
import importlib.util
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


tab = load("solid_brk_tab_check", os.path.join(HERE, "..", "solid-brk-tab", "check.py"))
RHO_S = 3660.0                                      # phase-file token rho-solid
A_V, B_V = 1.0 / tab.RHO_P, 1.0 / RHO_S - 1.0 / tab.RHO_P
M_0 = tab.RHO_P * math.pi / 6.0 * tab.D_0 ** 3      # both families inject liquid drops of d_0
HALF_ULP_F = 0.5e-6                                 # F10.6 f
DK, MK = "d<sub>p", "m<sub>p"


def hu_e13(v):
    """Half a unit in the last place of v printed E13.6 (mantissa 0.dddddd)."""
    return 0.5e-6 * 10.0 ** (math.floor(math.log10(abs(v))) + 1) if v else 0.0


def d_of(m, f):
    return (6.0 * m * (A_V + B_V * f) / math.pi) ** (1.0 / 3.0)


def main():
    try:
        traj, exits = tab.parcels(tab.TRAJ), tab.parcels(tab.EXIT)
        log = open(tab.LOG, errors="replace").read()
    except (FileNotFoundError, StopIteration, ValueError, KeyError) as e:
        print(f"[FAIL] cannot read the outputs: {e}")
        return 1
    print(f"oracle: m_0 = {M_0:.6e} kg, d_0 = {tab.D_0:.6e} m, rho_l = {tab.RHO_P}, rho_s = {RHO_S}: the frozen drop "
          f"is {(tab.RHO_P / RHO_S) ** (1.0 / 3.0):.6f} of the liquid one")
    fails = []

    # W1: the mass is fixed from the first frozen row on; family 2 keeps the injected one
    n1, w1_m0 = 0, 0.0
    for (fam, pid), rows in sorted(traj.items()):
        k = tab.first_frozen(rows)
        if k is None:
            fails.append(f"W1 family {fam} ID={pid}: no frozen row (phase >= 2) on its {len(rows)} rows")
            continue
        mk = rows[k][MK]
        bad = [r for r in rows[k + 1:] if r[MK] != mk]
        n1 += len(rows) - k
        if bad:
            fails.append(f"W1 family {fam} ID={pid}: {len(bad)} rows after the first frozen row (t={rows[k]['t']:.6e}, "
                         f"m={mk:.6e}) change m, first t={bad[0]['t']:.6e} m={bad[0][MK]:.6e} "
                         f"(ratio {bad[0][MK] / mk:.6f})")
        if fam == 2:
            dm = max(abs(r[MK] - M_0) for r in rows)
            w1_m0 = max(w1_m0, dm / hu_e13(M_0))
            if dm > hu_e13(M_0):
                fails.append(f"W1 family 2 ID={pid}: m departs from the injected {M_0:.6e} by {dm:.3e}")

    # W2: the diameter of the liquid and solid parts on every row
    w2, eff = (0.0, ""), {(fam, part): 0 for fam in (1, 2) for part in ("plateau", "solid")}
    n2 = 0
    for (fam, pid), rows in sorted(traj.items()):
        for r in rows:
            d, m, f = r[DK], r[MK], r["f"]
            tol = hu_e13(d) + d * (hu_e13(m) / m / 3.0 + abs(B_V / (3.0 * (A_V + B_V * f))) * HALF_ULP_F)
            res = abs(d - d_of(m, f))
            n2 += 1
            if res / tol > w2[0]:
                w2 = (res / tol, f"family {fam} ID={pid} t={r['t']:.6e} f={f:.6f} d={d:.6e} resid={res:.3e}")
            if res > tol:
                fails.append(f"W2 family {fam} ID={pid} t={r['t']:.6e}: d={d:.6e} vs (6 m (a + b f)/pi)^(1/3) = "
                             f"{d_of(m, f):.6e} at m={m:.6e}, f={f:.6f} (resid {res:.3e} > {tol:.3e})")
            elif f > 0.0 and abs(d - d_of(m, 0.0)) > tol:
                eff[(fam, "solid" if f >= 1.0 else "plateau")] += 1
    for key, n in eff.items():
        if n < tab.N_GOOD:
            fails.append(f"W2 non-vacuity: family {key[0]}: {n} {key[1]} rows where rho-solid moves d beyond the "
                         f"budget (need >= {tab.N_GOOD})")

    # W3: family 1 breaks up while liquid (by its mass)
    n3 = 0
    for (fam, pid), rows in sorted(traj.items()):
        if fam != 1:
            continue
        k = tab.first_frozen(rows)
        dec = [i for i in range(1, len(rows)) if rows[i][MK] < rows[i - 1][MK]]
        if k is None or not any(i <= k for i in dec):
            fails.append(f"W3 family 1 ID={pid}: no m decrease up to its first frozen row (k={k}, decreases at {dec[:4]})")
            continue
        n3 += 1

    w4 = tab.check_flows(exits, fails, "T4", per_record=True)
    got, ref, rel5 = tab.check_energy(exits, fails, "T5")
    good = tab.check_injection(traj, exits, fails, "T6")
    if log.count("solidification --> supercooling + recalescence (model 6)") != 2:
        fails.append("T6 witness: the log does not report the model-6 ODE system for both groups")

    print(f"W1 {n1} rows from the first frozen row on; family 2 worst |m - m_0| = {w1_m0:.3f} half-ULP")
    print(f"W2 {n2} rows: worst d resid/tol {w2[0]:.3f} ({w2[1]}); rows where rho-solid moves d: "
          + ", ".join(f"family {k[0]} {k[1]} {n}" for k, n in sorted(eff.items())))
    print(f"W3 {n3} family-1 parcels lose mass by breakup while liquid")
    print(f"T4 {sum(len(v) for v in exits.values())} exit records; worst exit-flow deviation {w4:.3e} (tol {tab.TOL_E:.0e})")
    print(f"T5 sum(E) = {got:.9e} W, closed form {ref:.9e} W, rel {rel5:.3e}")
    print(f"T6 parcels verified per family: {good}")
    for f in fails[:40]:
        print(f"[FAIL] {f}")
    if fails:
        print(f"\n[FAIL] {len(fails)} violation(s)")
        return 1
    print("\n[PASS] the frozen drop keeps the mass of its last breakup and takes the diameter of its liquid and solid "
          "parts; flows and energy balance.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
