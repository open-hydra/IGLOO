#!/usr/bin/env python3
"""The check shared by body-evap-m2, body-evap-m4 and body-burn-m5: the body-force accumulators of parcels whose mass
changes (ODE models 2, 4 and 5), seen through the coupling source.

A body acceleration g gives a parcel momentum and energy the gas never gave it. The source hands the gas the
reaction of the drag and the mass the parcel loses, so each segment's deposit carries the correction
P += ndot J g and E += ndot W, J = int m dt and W = int m g.v dt integrated along the path as extra ODE states
(models 2 and 5: J and W at neq-1 and neq; model 4, whose states carry the stream's number rate: J at neq-2 when the
Eulerian output is on, neq-1 then holding the number rate). Momentum and energy telescope along each path:

  B1  sum Fy = sum_p ndot [(m v_y)_first - (m v_y)_last] + g_y sum_p ndot int m dt
  B2  sum E  = sum_p ndot [(m (h + |v|^2/2))_first - (m (h + |v|^2/2))_last] + g_y sum_p ndot int m v_y dt

read from the trajectory rows, which the source path does not write: first = the injection state, last = the state
at time-end (every parcel stops mid-box before any reaches a wall, so its last row is its final state), ndot the
injected number rate krho/(1 - krho) rho_g U dy dz / m_0, h = c_p T (a constant-cp table, relative datum). The
integrals are trapezoids over the printed time (nine digits), each with a Richardson estimate (all rows against
every other row). The enthalpy flux m c_p T is about 1e6 times the body work and the E13.6 print of m alone is
about 1e-6 of it, so B2 values the mass the parcels lost by the mass source itself: sum_p ndot m_last h_last =
h_last (sum_p ndot m_0 - sum wdot) -- the 25 parcels ride one uniform gas from one injection state and end in one
state within the integrator's spread (asserted; h_last at their mean T, the spread in B2's tolerance), and B0 holds
sum wdot to the rows.

Gates, on the run in the case directory (out-file = ALL) and on a sibling ref-S/ with out-file = S (the two slot
layouts of the accumulators):
  W   witnesses: the body acceleration, the output fields, the ODE size of the layout; 25 parcels, each with one exit
      record, code 0 (none left through a face) at its snapshot row's position and time; each first row the injection
      state and each last row at t = time-end exactly, the snapshot holding them; the
      parcels' final T within SPREAD_T of their mean and m within SPREAD_REL (B2's premise only)
  B0  sum wdot = sum_p ndot (m_0 - m_last) and sum Fx = sum_p ndot u_g (m_0 - m_last) (no body term), within the
      print of m; sum Fz = 0
  B1  within 1e-5 of the gross plus 3 Richardson estimates
  B2  within the prints of T and of the kinetic term, the final states' spread and 3 Richardson estimates
  B3  (evaporation) d^2 along the d2-law kernel 8 k_g/(c_p,g rho_p) ln(1 + c_p,g (T_g - T)/L_v) integrated over the
      measured T (the stagnant-film law does not depend on the slip), within the prints and the Richardson estimate
  B4  non-vacuity: g_y sum_p ndot int (m_0 - m) dt, what a closed form J = m_0 t would add, exceeds 100 B1 tolerances
  D   (model 4) the breakup stays dormant: We < 6 and We/Re^1/2 < 0.5 on every row, and sum n_p V of the Eulerian
      field equals 25 ndot time-end (the number rate never moved)
  B5  (a regression band, no oracle) the two runs' source.tec agree cell by cell within BAND of each variable's
      scale: the two layouts integrate one ODE under different error norms
With IGLOO_REUSE_SIBLING=1 the ref-S/ a previous run left is read, not run again.
"""
import math
import os
import re
import shutil
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import datcols  # noqa: E402

# ---- the box (make_box_case.py defaults) and the shared inputs (SI) --------------------------------------------
RHO_G, U_G, T_G, MU_G, K_G, GAM, R_G = 1.2, 10.0, 600.0, 1.8e-5, 0.026, 1.4, 287.0
CP_G = GAM * R_G / (GAM - 1.0)
KRHO, DY, DZ, CELL_V = 0.34, 0.01, 0.01, 0.0025 * 0.01 * 0.01   # bc.txt krho; inlet cell; box cell volume
G_Y, T_END, N_PARCELS = -200.0, 6.875e-3, 25                       # input.ini body-accel, time-end

# ---- error model ---------------------------------------------------------------------------------------------
HALF_F = 0.5e-6           # F12.6 half-unit (X, U, V, T)
TOL_B1_REL = 1.0e-5       # of the gross, plus 3 Richardson estimates
INT_FLOOR_REL = 1.0e-7    # B3: the integrator, relative to d_0^2
SUM_FLOOR = 1.0e-12       # relative: the 15-digit print of the source fields and the sums
SPREAD_T, SPREAD_REL = 5.0e-5, 1.0e-5   # B2's premise: the parcels' final T [K] and m (relative) about their mean
SIGNAL = 100.0            # B4

EXITS, TRAJ, SNAP, SRC, EUL, LOG = ("OUTPUT/outloc-A.dat", "OUTPUT/trajectories-A.dat", "OUTPUT/snapshot-A.dat",
                                     "OUTPUT/source.tec", "OUTPUT/euler1.tec", "run_out.txt")


def half_e(v):
    """Half a unit of the last digit of an E13.6 print (0.dddddd E+xx)."""
    return 0.5e-6 * 10.0 ** (math.floor(math.log10(abs(v))) + 1) if v else 0.0


def half_es(t):
    """Half a unit of the last digit of an ES16.8 print (d.dddddddd E+xx)."""
    return 0.5e-8 * 10.0 ** math.floor(math.log10(abs(t))) if t else 0.0


def cells(path):
    """The cell-centred variables of a BLOCK-packed Tecplot file, by the names of its VARIABLES line."""
    lines = open(path).read().splitlines()
    names = lines[0].split("=", 1)[1].split()
    zi = next(n for n, l in enumerate(lines) if l.strip().upper().startswith("ZONE"))
    I, J, K = (int(re.search(r"\b%s=(\d+)" % a, lines[zi]).group(1)) for a in "IJK")
    nn, nc = I * J * K, (I - 1) * (J - 1) * (K - 1)
    vals = [float(v) for v in " ".join(lines[zi + 1:]).split()]
    return {names[3 + a]: vals[3 * nn + a * nc: 3 * nn + (a + 1) * nc] for a in range(len(names) - 3)}


def trap(ts, fs):
    return sum(0.5 * (fs[i] + fs[i + 1]) * (ts[i + 1] - ts[i]) for i in range(len(ts) - 1))


def trap_rich(ts, fs):
    """(trapezoid over all rows, |all rows - every other row|): the coarse grid keeps the last row."""
    idx = list(range(0, len(ts) - 1, 2)) + [len(ts) - 1]
    fine = trap(ts, fs)
    return fine, abs(fine - trap([ts[i] for i in idx], [fs[i] for i in idx]))


def run_sibling(here, igloo):
    """ref-S/: the case with out-file = S (the source alone), INPUT linked."""
    ref = os.path.join(here, "ref-S")
    shutil.rmtree(ref, ignore_errors=True)
    os.makedirs(os.path.join(ref, "OUTPUT"))
    os.symlink(os.path.join(here, "INPUT"), os.path.join(ref, "INPUT"))
    text = open(os.path.join(here, "input.ini")).read()
    new, n = re.subn(r"^out-file\s*=\s*ALL\s*$", "out-file    = S", text, flags=re.M)
    if n != 1:
        raise RuntimeError("input.ini: no 'out-file = ALL' line to switch")
    open(os.path.join(ref, "input.ini"), "w").write(new)
    with open(os.path.join(ref, "run_out.txt"), "w") as out, open(os.path.join(ref, "run_err.txt"), "w") as err:
        return subprocess.run([igloo], cwd=ref, stdout=out, stderr=err).returncode


class Gates:
    def __init__(self, tag):
        self.tag, self.fails, self.lines = tag, [], []

    def check(self, gate, ok, msg):
        line = "[%s] %s %s: %s" % ("PASS" if ok else "FAIL", self.tag, gate, msg)
        self.lines.append(line)
        if not ok:
            self.fails.append(line)


def check_run(d, cfg, layout, g):
    """Every gate of one run in directory d; returns the source fields for B5."""
    m0 = cfg["rho_p"] * math.pi / 6.0 * cfg["d0"] ** 3
    ndot = KRHO / (1.0 - KRHO) * RHO_G * U_G * DY * DZ / m0
    log = open(os.path.join(d, LOG), errors="replace").read()
    field = ("Output fields: gas coupling source & equivalent eulerian" if layout == "ALL"
             else "Output field: gas coupling source")
    neq = cfg["neq_" + layout]
    g.check("W", "Body acceleration ON" in log and field in log and "ODE system (neq=%d)" % neq in log,
            "log: body acceleration, '%s', ODE size %d" % (field, neq))
    exits = datcols.rows(os.path.join(d, EXITS))
    rows = datcols.rows(os.path.join(d, TRAJ))
    snap = {r["ID"]: r for r in datcols.rows(os.path.join(d, SNAP))}
    src = cells(os.path.join(d, SRC))
    parts = {}
    for r in rows:
        parts.setdefault(r["ID"], []).append(r)
    for p in parts.values():
        p.sort(key=lambda r: r["t"])
    parts = dict(sorted(parts.items()))  # ID order: the sums and the ties do not follow the record order
    ex_ok = sorted(r["ID"] for r in exits) == sorted(snap) and all(
        r["bc"] == 0 and all(r[k] == snap[r["ID"]][k] for k in ("X", "Y", "Z", "t")) for r in exits)
    g.check("W", ex_ok and len(exits) == N_PARCELS and len(parts) == N_PARCELS and len(snap) == N_PARCELS,
            "%d exit records, each code 0 at its snapshot row's position and time (time-end comes first), %d parcels, "
            "%d snapshot rows" % (len(exits), len(parts), len(snap)))
    keys = ("X", "U", "V", "W", "T", "d<sub>p", "m<sub>p", "t")
    bad_first = [pid for pid, p in parts.items()
                 if p[0]["t"] != 0.0 or abs(p[0]["T"] - cfg["T0"]) > HALF_F or abs(p[0]["U"] - U_G) > HALF_F
                 or abs(p[0]["V"]) > HALF_F or abs(p[0]["m<sub>p"] - m0) > half_e(m0)]
    bad_last = [pid for pid, p in parts.items()
                if p[-1]["t"] != T_END or pid not in snap or any(snap[pid][k] != p[-1][k] for k in keys)]
    g.check("W", not bad_first and not bad_last,
            "first rows the injection state (T0 %.1f K, u_g, m_0 %.6e kg), last rows at t = %.6e s = the snapshot"
            % (cfg["T0"], m0, T_END))
    if bad_first or bad_last:
        return src
    t_bar = sum(p[-1]["T"] for p in parts.values()) / len(parts)
    m_bar = sum(p[-1]["m<sub>p"] for p in parts.values()) / len(parts)
    spread_t = max(abs(p[-1]["T"] - t_bar) for p in parts.values())
    spread_m = max(abs(p[-1]["m<sub>p"] / m_bar - 1.0) for p in parts.values())
    one = spread_t <= SPREAD_T and spread_m <= SPREAD_REL
    g.check("W", one, "one final state (B2's premise): the parcels' final T within %.1e K of their mean, m within %.1e"
            % (spread_t, spread_m))

    # ---- B0: mass and x-momentum, no body term ----
    lost = sum(ndot * (m0 - p[-1]["m<sub>p"]) for p in parts.values())
    tol0 = sum(ndot * half_e(p[-1]["m<sub>p"]) for p in parts.values()) + SUM_FLOOR * ndot * m0 * N_PARCELS
    wsum, fxsum, fzsum = sum(src["wdot(A)"]), sum(src["Fx"]), sum(src["Fz"])
    g.check("B0", abs(wsum - lost) <= tol0, "sum wdot %.9e vs rows %.9e kg/s (%.3f of the m print)" % (
        wsum, lost, abs(wsum - lost) / tol0))
    g.check("B0", abs(fxsum - U_G * lost) <= U_G * tol0 + HALF_F * lost,
            "sum Fx %.9e vs u_g x the mass lost %.9e N (%.3f)" % (
                fxsum, U_G * lost, abs(fxsum - U_G * lost) / (U_G * tol0 + HALF_F * lost)))
    g.check("B0", abs(fzsum) <= SUM_FLOOR * abs(fxsum), "sum Fz %.3e (no motion along z)" % fzsum)

    # ---- B1 / B2 / B4 ----
    b1 = gross1 = ric1 = 0.0
    b2_first = b2_kin = wint = ric2 = tol_t2 = tol_k2 = 0.0
    sig4 = 0.0
    for p in parts.values():
        ts = [r["t"] for r in p]
        ms = [r["m<sub>p"] for r in p]
        im, rm = trap_rich(ts, ms)
        imv, rmv = trap_rich(ts, [r["m<sub>p"] * r["V"] for r in p])
        last = p[-1]
        b1 += ndot * (0.0 - last["m<sub>p"] * last["V"]) + G_Y * ndot * im
        gross1 += ndot * (abs(last["m<sub>p"] * last["V"]) + abs(G_Y * im))
        ric1 += abs(G_Y * ndot * rm)
        v2 = last["U"] ** 2 + last["V"] ** 2 + last["W"] ** 2
        b2_first += ndot * m0 * (cfg["cp"] * cfg["T0"] + 0.5 * U_G * U_G)
        b2_kin += ndot * last["m<sub>p"] * 0.5 * v2
        wint += G_Y * ndot * imv
        ric2 += abs(G_Y * ndot * rmv)
        tol_t2 += ndot * last["m<sub>p"] * cfg["cp"] * HALF_F
        tol_k2 += ndot * (half_e(last["m<sub>p"]) * 0.5 * v2 + last["m<sub>p"] * math.sqrt(v2) * 2.0 * HALF_F)
        sig4 += G_Y * ndot * (m0 * T_END - im)
    tol1 = TOL_B1_REL * gross1 + 3.0 * ric1
    fysum = sum(src["Fy"])
    g.check("B1", abs(fysum - b1) <= tol1, "sum Fy %.9e vs telescoping %.9e N: %.3e (%.3f of tol %.3e)" % (
        fysum, b1, abs(fysum - b1), abs(fysum - b1) / tol1, tol1))
    h_last = cfg["cp"] * t_bar
    b2 = b2_first - h_last * (N_PARCELS * ndot * m0 - wsum) - b2_kin + wint
    tol2 = (tol_t2 * (1.0 + spread_t / HALF_F) + tol_k2 + 3.0 * ric2 + SUM_FLOOR * b2_first)
    esum = sum(src["E"])
    g.check("B2", one and abs(esum - b2) <= tol2, "sum E %.9e vs telescoping %.9e W: %.3e (%.3f of tol %.3e; the body "
            "work %.3e W is %.0f tolerances)" % (esum, b2, abs(esum - b2), abs(esum - b2) / tol2, tol2, wint,
                                                 abs(wint) / tol2))
    g.check("B4", abs(sig4) > SIGNAL * tol1, "g_y sum ndot int (m_0 - m) dt = %.3e N, %.0f B1 tolerances" % (
        sig4, abs(sig4) / tol1))

    # ---- B3: the d2-law kernel along the measured temperature ----
    if cfg.get("lv"):
        apre = 8.0 * K_G / (CP_G * cfg["rho_p"])
        worst, nbad = (0.0, ""), 0
        for pid, p in parts.items():
            ts = [r["t"] for r in p]
            ks = [apre * math.log1p(max(CP_G * (T_G - r["T"]) / cfg["lv"], 0.0)) for r in p]
            dks = [apre * CP_G / (cfg["lv"] + CP_G * (T_G - r["T"])) for r in p]
            d2a = p[0]["d<sub>p"] ** 2
            for k in range(1, len(p)):
                integ, rich = trap_rich(ts[:k + 1], ks[:k + 1])
                eT = trap(ts[:k + 1], dks[:k + 1]) * HALF_F
                dk = p[k]["d<sub>p"]
                tol = (2.0 * dk * half_e(dk) + 2.0 * p[0]["d<sub>p"] * half_e(p[0]["d<sub>p"]) + rich + eT
                       + INT_FLOOR_REL * d2a)
                res = abs(dk * dk - (d2a - integ))
                if res / tol > worst[0]:
                    worst = (res / tol, "ID=%d t=%.6e" % (pid, ts[k]))
                nbad += res > tol
        g.check("B3", nbad == 0, "d^2 on the d2-law kernel along the measured T: %d row(s) off, worst %.3f of "
                "the tolerance (%s)" % (nbad, worst[0], worst[1]))

    # ---- D: model 4's breakup stays dormant ----
    if cfg.get("sigma"):
        wmax = rmax = 0.0
        for p in parts.values():
            for r in p:
                slip = math.hypot(r["U"] - U_G, r["V"])
                we = RHO_G * slip * slip * r["d<sub>p"] / cfg["sigma"]
                re_ = RHO_G * slip * r["d<sub>p"] / MU_G
                wmax = max(wmax, we)
                rmax = max(rmax, we / math.sqrt(re_) if re_ > 0.0 else 0.0)
        g.check("D", wmax < 6.0 and rmax < 0.5, "We max %.3e (< 6), We/Re^1/2 max %.3e (< 0.5)" % (wmax, rmax))
        if layout == "ALL":
            ntot = sum(cells(os.path.join(d, EUL))["n<sub>p"]) * CELL_V
            nref = N_PARCELS * ndot * T_END
            g.check("D", abs(ntot - nref) <= SUM_FLOOR * nref * 1e3,
                    "sum n_p V %.12e vs 25 ndot time-end %.12e (the number rate never moved)" % (ntot, nref))
    return src


def main(cfg):
    here = os.getcwd()
    if os.environ.get("IGLOO_REUSE_SIBLING") == "1":  # read the ref-S/ a previous run left
        if not os.path.isdir(os.path.join(here, "ref-S", "OUTPUT")):
            print("[FAIL] IGLOO_REUSE_SIBLING=1 and no ref-S/OUTPUT to read")
            return 1
    else:
        igloo = os.environ.get("IGLOO_EXE", os.path.abspath(os.path.join(here, "..", "..", "..", "bin", "IGLOO")))
        rc = run_sibling(here, igloo)
        if rc != 0:
            print("[FAIL] the ref-S sibling run failed (exit %d)" % rc)
            return 1
    ga, gs = Gates("ALL"), Gates("ref-S")
    try:
        sa = check_run(here, cfg, "ALL", ga)
        ss = check_run(os.path.join(here, "ref-S"), cfg, "S", gs)
    except (OSError, ValueError, KeyError, StopIteration) as exc:
        print("[FAIL] cannot read the outputs: %s" % exc)
        return 1
    gb = Gates("pair")
    worst = (0.0, "")
    for k in ("wdot(A)", "Fx", "Fy", "E"):
        scale = max(abs(v) for v in sa[k])
        diff = max(abs(a - b) for a, b in zip(sa[k], ss[k]))
        if diff / scale > worst[0]:
            worst = (diff / scale, k)
    gb.check("B5", worst[0] <= cfg["band"], "ALL and S source.tec cell by cell: worst %.3e of the scale (%s), "
             "band %.1e (a regression band)" % (worst[0], worst[1], cfg["band"]))
    for g in (ga, gs, gb):
        print("\n".join(g.lines))
    fails = ga.fails + gs.fails + gb.fails
    if fails:
        print("\n[FAIL] %d gate(s) failed" % len(fails))
        return 1
    print("\n[PASS] model %d: the body-force accumulators close momentum and energy in both slot layouts" % cfg["model"])
    return 0
