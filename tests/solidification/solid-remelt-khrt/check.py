#!/usr/bin/env python3
"""solid-remelt-khrt: a KH-RT drop that freezes, re-melts and breaks up (ODE model 7).

solid-brk-khrt's material and drop (30 um, 35 m/s) enter the box gas (10 m/s) 1 K above T_nuc, with the gas
temperature rising along x (T = 1500 + 20000 x K at the cell centres). The drop nucleates within a few microseconds,
its plateau first freezes, then melts once the gas is hotter than T_m, and re-melts (f = 0) near the end of the box.
The re-melt restarts the KH-RT clock: the first RT breakup follows the re-melt by tau_RT, the Beale-Reitz time of the
drop at the event. KH stripping needs We_r above WeLimit (6); at the re-melt the drop's We_r is 4.0, so its d stays
the re-melt row's d up to the RT.

Rows are read in time order per parcel (tools/datcols.py); a row's phase is datcols.solid_phase(f, T, T_m).

R0  per parcel: a row in phase 2 or 3 before 0.1 ms; a later row with f = 0 in phase 0 (the re-melt, t_rm); no row
    before t_rm with d below 0.9 d_0; the first row after t_rm with d below 0.9 of the previous row's (the RT, t_RT),
    every row between them with the re-melt row's d; lambda_RT < d at both rows (the clock runs throughout); the RT
    row's d is lambda_RT of the re-melt row's d at the RT row's u_p (one RT breakup, its products of diameter
    lambda_RT), within the E13.6 half-ULPs of the two d. Every injected parcel (one per inlet record of bc.txt) is
    checked. Exact on the printed values but the product diameter.
R1  per parcel: |(t_RT - t_rm) - tau_RT| <= 1e-8 s, tau_RT = Ctau/Omega_RT from the RT row's u_p and the re-melt
    row's d: a = 18 mu_g (u_g - u_p)/(rho_l d^2) (the case's Stokes law), its component along u_p, force =
    |a_t| (rho_l - rho_g), Omega_RT = (2 force^(3/2)/(3 (3 sigma)^(1/2) (rho_l + rho_g)))^(1/2).
    Inputs from input.ini (sigma, Ctau, CRT), properties.dat (Density), solfile.tec (the uniform gas density,
    velocity and viscosity), phase.txt (T-melt), bc.txt (the injected radius).
"""
import configparser
import math
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "tools"))
import datcols                                  # noqa: E402

TRAJ = "OUTPUT/trajectories-A.dat"
T_NUC_BY = 1.0e-4         # R0: the nucleation row before 0.1 ms
D_DROP = 0.9              # R0: a breakup row has d below 0.9 of the previous row's
TOL_R1 = 1.0e-8           # R1 [s]


def ini():
    cp = configparser.ConfigParser(comment_prefixes=(";", "#"), inline_comment_prefixes=(";",), strict=False)
    cp.optionxform = str
    cp.read(os.path.join(HERE, "input.ini"))
    models = cp["IGLOO-Models"]
    return (float(cp["IGLOO-Properties"]["sigma"]), float(models.get("Ctau", "1")),
            float(models.get("CRT", "0.1")))


def table_column(path, name):
    """The values of one column of a Tecplot POINT table."""
    lines = open(path).read().splitlines()
    cols = datcols._names(next(l for l in lines if l.strip().lower().startswith("variables")))
    k = cols.index(name)
    out = []
    for l in lines:
        v = l.split()
        if len(v) == len(cols):
            try:
                out.append(float(v[k]))
            except ValueError:
                pass
    return out


def gas():
    """(rho_g, u_g, mu_g) of the uniform box gas, from the cell-centred block of solfile.tec."""
    lines = open(os.path.join(HERE, "INPUT", "solfile.tec")).read().splitlines()
    names = re.findall(r'"([^"]*)"', next(l for l in lines if "VARIABLES" in l.upper()))
    zi = next(i for i, l in enumerate(lines) if l.strip().upper().startswith("ZONE"))
    I, J, K = (int(re.search(r"\b%s=(\d+)" % a, lines[zi]).group(1)) for a in "IJK")
    vals = [float(t) for l in lines[zi + 2:] for t in l.split()]
    nn, nc = I * J * K, (I - 1) * (J - 1) * (K - 1)
    out = {}
    for n in ("Roi( 1 )", "U", "V", "W", "MIL"):
        k = names.index(n) - 3
        blk = vals[3 * nn + k * nc: 3 * nn + (k + 1) * nc]
        if min(blk) != max(blk):
            raise SystemExit("[FAIL] solfile.tec: %s is not uniform" % n)
        out[n] = blk[0]
    if out["V"] != 0.0 or out["W"] != 0.0:
        raise SystemExit("[FAIL] solfile.tec: the gas moves off x")
    return out["Roi( 1 )"], out["U"], out["MIL"]


def phase_tmelt():
    toks = open(os.path.join(HERE, "INPUT", "phase.txt")).read().split()
    tm = [float(t.split("=")[1]) for t in toks if t.lower().startswith("t-melt=")]
    return tm[0] if tm else 2327.0          # the T-melt default


def inlet():
    """(number of inlet records, the injected diameter) from bc.txt: a 401 header, then its property line."""
    lines = open(os.path.join(HERE, "INPUT", "bc.txt")).read().splitlines()
    n, rp = 0, set()
    for i, l in enumerate(lines):
        v = l.split()
        if len(v) == 6 and v[5] == "401":
            n += 1
            rp.add(float(lines[i + 1].split()[5]))
    if len(rp) != 1:
        raise SystemExit("[FAIL] bc.txt: %d inlet radii" % len(rp))
    return n, 2.0 * rp.pop()


def hulp_e13(v):
    """Half a unit in the last place of an E13.6 print (0.dddddd mantissa) of v > 0."""
    return 0.5e-6 * 10.0 ** (math.floor(math.log10(v)) + 1)


def waves(d, up, rho_l, rho_g, u_g, mu_g, sigma, ctau, crt):
    """(tau_RT, lambda_RT) of a drop of diameter d at velocity up (a 3-vector) in the uniform gas."""
    k = 18.0 * mu_g / (rho_l * d * d)
    acc = [k * (u_g - up[0]), -k * up[1], -k * up[2]]
    at = sum(a * u for a, u in zip(acc, up)) / math.sqrt(sum(u * u for u in up))
    force = abs(at) * (rho_l - rho_g)
    omega = math.sqrt(2.0 * force ** 1.5 / (3.0 * math.sqrt(3.0 * sigma) * (rho_l + rho_g)))
    return ctau / omega, 2.0 * math.pi * crt / math.sqrt(force / (3.0 * sigma))


def main():
    os.chdir(HERE)
    sigma, ctau, crt = ini()
    rho = table_column(os.path.join("INPUT", "properties.dat"), "Density")
    if min(rho) != max(rho):
        raise SystemExit("[FAIL] properties.dat: the density varies (this oracle takes it constant)")
    rho_l = rho[0]
    rho_g, u_g, mu_g = gas()
    t_m = phase_tmelt()
    n_inj, d_0 = inlet()

    parcels = {}
    for z in datcols.parcel_zones(TRAJ):
        parcels.setdefault((z["group"], z["id"]), []).extend(z["rows"])
    fails = []
    if len(parcels) != n_inj:
        fails.append("R0 %d parcels in the trajectories, %d inlet records" % (len(parcels), n_inj))
    worst, prod, lam, seq = 0.0, 0.0, [], []
    for key, rows in sorted(parcels.items()):
        rows.sort(key=lambda r: r["t"])
        ph = [datcols.solid_phase(r["f"], r["T"], t_m) for r in rows]
        tag = "parcel %d/%d" % key
        k_n = next((i for i, r in enumerate(rows) if ph[i] >= 2), None)
        if k_n is None or rows[k_n]["t"] >= T_NUC_BY:
            fails.append("R0 %s: no frozen row before %.1e s" % (tag, T_NUC_BY))
            continue
        k_m = next((i for i in range(k_n + 1, len(rows)) if rows[i]["f"] == 0.0 and ph[i] == 0), None)
        if k_m is None:
            fails.append("R0 %s: no re-melt row" % tag)
            continue
        low = [i for i in range(k_m) if rows[i]["d<sub>p"] < D_DROP * d_0]
        if low:
            fails.append("R0 %s: row %d before the re-melt with d %.6e" % (tag, low[0], rows[low[0]]["d<sub>p"]))
        d_m = rows[k_m]["d<sub>p"]
        k_r = next((i for i in range(k_m + 1, len(rows)) if rows[i]["d<sub>p"] < D_DROP * rows[i - 1]["d<sub>p"]),
                   None)
        if k_r is None:
            fails.append("R0 %s: no breakup row after the re-melt" % tag)
            continue
        moved = [i for i in range(k_m + 1, k_r) if rows[i]["d<sub>p"] != d_m]
        if moved:
            fails.append("R0 %s: row %d between the re-melt and the RT with d %.6e (re-melt d %.6e)"
                         % (tag, moved[0], rows[moved[0]]["d<sub>p"], d_m))
        args = (rho_l, rho_g, u_g, mu_g, sigma, ctau, crt)
        _, lam_m = waves(d_m, [rows[k_m][c] for c in ("U", "V", "W")], *args)
        tau, lam_r = waves(d_m, [rows[k_r][c] for c in ("U", "V", "W")], *args)
        if not (lam_m < d_m and lam_r < d_m):
            fails.append("R0 %s: lambda_RT/d %.4f, %.4f (the clock stops)" % (tag, lam_m / d_m, lam_r / d_m))
        lam.append(max(lam_m, lam_r) / d_m)
        d_r = rows[k_r]["d<sub>p"]
        bound = hulp_e13(d_r) + lam_r / d_m * hulp_e13(d_m)
        prod = max(prod, abs(d_r - lam_r) / bound)
        if abs(d_r - lam_r) > bound:
            fails.append("R0 %s: the RT row's d %.6e, lambda_RT %.6e (off by %.2e > %.2e: not one RT breakup)"
                         % (tag, d_r, lam_r, abs(d_r - lam_r), bound))
        res = abs((rows[k_r]["t"] - rows[k_m]["t"]) - tau)
        worst = max(worst, res)
        if res > TOL_R1:
            fails.append("R1 %s: t_RT - t_rm = %.6e s, tau_RT = %.6e s, off by %.2e > %.0e"
                         % (tag, rows[k_r]["t"] - rows[k_m]["t"], tau, res, TOL_R1))
        seq.append((rows[k_n]["t"], rows[k_m]["t"], rows[k_m]["X"], rows[k_r]["t"], rows[k_r]["X"], tau))
    if seq:
        tn, tm_, xm, tr, xr, tau = seq[0]
        print("[INFO] first parcel: nucleation %.4e s; re-melt %.6e s at x %.6f; RT %.6e s at x %.6f; tau_RT %.6e s"
              % (tn, tm_, xm, tr, xr, tau))
        print("[INFO] R0: %d parcels freeze, re-melt and break up; lambda_RT/d at most %.4f; the RT row's d against "
              "lambda_RT at worst %.3f of its bound" % (len(seq), max(lam), prod))
        print("[INFO] R1: worst |(t_RT - t_rm) - tau_RT| = %.3e s against %.0e" % (worst, TOL_R1))
    for f in fails:
        print("[FAIL] " + f)
    if fails or not seq:
        print("solid-remelt-khrt: FAIL")
        return 1
    print("solid-remelt-khrt: PASS (R0 R1)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
