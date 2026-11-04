#!/usr/bin/env python3
"""tab-hf92 (etab-hf92 runs this file on its own run): TAB or ETAB products of Hsiang & Faeth's water series.

Every inlet cell injects d0 = 1000 um water at its own slip u0, one diameter-based Weber number per cell
(We = rho_G d0 u0^2/sigma, 15 to 600): a class is an inlet cell. verify.py draws the classes against [HF92]
eq. [14] (the SMD), the root-normal size distribution ([HF92] p. 645, [HF93] Fig. 1) and [HF93] eq. [13]
(drop size against relative velocity at the end of breakup, t_b = 5 t*); those comparisons are drawn, never
gated. This oracle gates what the panels need:
  P0  every class holds >= 15 exit records; a parcel's class is the inlet cell holding the (Y, Z) of its
      injection row (t = 0), whose |u_p| must be that cell's kv U; one exit record per parcel
  P1  every parcel leaves with d < 0.999 d0 (every class is above the bag onset We = 13)
  P2  every parcel leaves through the outlet (X >= Lx - 1e-6), keeps one d over the box's last half, and the
      radius Weber number of its last row, rho_G |u_g - u_p|^2 (d/2)/sigma, is below 6 (TAB and ETAB break
      only above it)
  P3  every parcel's rows bracket t_b = 5 t* of its class, t* = d0 (rho_L/rho_G)^1/2/u0
Inputs come from the case files: kv and d0 from INPUT/bc.txt, the gas and the box from INPUT/solfile.tec,
sigma and mu from input.ini, the density from INPUT/properties.dat. Parcels are taken by ascending ID and a
parcel's rows by time, so no result depends on the order of the records.
"""
import math
import os
import re
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "tools"))
from datcols import rows  # noqa: E402
from proptab import load_properties, table_value  # noqa: E402

N_MIN = 15          # P0: exit records per class
BROKE = 0.999       # P1: d_exit/d0 below it
WE_R_MAX = 6.0      # P2: radius Weber number of the last row below it
X_TOL = 1.0e-6      # P2: the outlet plane at the F12.6 print resolution; P0: the cell edges likewise
U_TOL = 1.0e-5      # P0: |u_p| of the injection row against kv U (F12.6)
TB_STAR = 5.0       # t_b/t* ([HF92] eqs. [1]-[2], [HF93] p. 726)
C_SMD = 6.2         # [HF92] eq. [14]
ONSETS = ((13.0, "bag"), (35.0, "multimode"), (80.0, "shear"))   # [HF93] p. 724
INT = re.compile(r"[+-]?\d+$")


def read_ini(path):
    sec, out = None, {}
    for ln in open(path):
        s = ln.strip()
        if not s or s[0] in ";#!":
            continue
        if s.startswith("[") and s.endswith("]"):
            sec = s[1:-1]
            continue
        if "=" in s and sec is not None:
            k, v = s.split("=", 1)
            out[(sec, k.strip())] = re.split(r"[;#!]", v)[0].strip()
    return out


def read_solfile(path):
    """Box length, inlet-cell edges in Y and Z, and the uniform gas of a make_pe_case.py solfile.tec."""
    with open(path) as fh:
        head = [fh.readline() for _ in range(4)]
        vals = fh.read().split()
    names = re.findall(r'"([^"]*)"', head[1])
    ni, nj, nk = (int(re.search(r"\b%s=(\d+)" % c, head[2]).group(1)) for c in "IJK")
    nn, nc = ni * nj * nk, (ni - 1) * (nj - 1) * (nk - 1)
    lx = max(float(v) for v in vals[:nn])
    ys = [float(vals[nn + ni * j]) for j in range(nj)]
    zs = [float(vals[2 * nn + ni * nj * k]) for k in range(nk)]
    gas = {}
    for i, n in enumerate(names[3:]):
        blk = vals[3 * nn + i * nc:3 * nn + (i + 1) * nc]
        if len(set(blk)) != 1:
            raise SystemExit("[FAIL] %s: %s is not uniform; the classes assume one gas state" % (path, n))
        gas[n] = float(blk[0])
    return lx, ys, zs, gas


def read_bc(path):
    """[(m, n, kv, rp)] of the face-1 inlet (401) records, in file order."""
    lines = open(path).read().splitlines()
    out = []
    for a, b in zip(lines, lines[1:]):
        h = a.split()
        if len(h) == 6 and all(INT.match(t) for t in h) and int(h[1]) == 1 and int(h[5]) == 401:
            p = b.split()
            out.append((int(h[2]), int(h[3]), float(p[1]), float(p[5])))
    return out


def dkey(row):
    return next(k for k in row if k.startswith("d<sub>p"))


def speed(r):
    return math.sqrt(r["U"] ** 2 + r["V"] ** 2 + r["W"] ** 2)


def hf92_smd(c, case):
    """[HF92] eq. [14]: SMD/d0 = 6.2 (rho_L/rho_G)^1/4 [mu_L/(rho_L d0 u0)]^1/2."""
    return C_SMD * (case["rhol"] / case["rhog"]) ** 0.25 * math.sqrt(case["mul"] / (case["rhol"] * c["d0"] * c["u0"]))


def mmd(ds, ws):
    """The d at cumulative volume fraction 0.5, linear between the parcels' mid-step positions."""
    pts = sorted(zip(ds, ws))
    tot, acc, pos = sum(w for _, w in pts), 0.0, []
    for d, w in pts:
        pos.append((d, (acc + 0.5 * w) / tot))
        acc += w
    if pos[0][1] >= 0.5:
        return pos[0][0]
    for (d0, v0), (d1, v1) in zip(pos, pos[1:]):
        if v1 >= 0.5:
            return d0 + (d1 - d0) * (0.5 - v0) / (v1 - v0) if v1 > v0 else d1
    return pos[-1][0]


def at_time(seq, t):
    """X, d and the slip vector at time t, linear between the rows bracketing it (None outside them)."""
    for a, b in zip(seq, seq[1:]):
        if a["t"] <= t <= b["t"] and b["t"] > a["t"]:
            f = (t - a["t"]) / (b["t"] - a["t"])
            g = {k: a[k] + f * (b[k] - a[k]) for k in ("X", "U", "V", "W")}
            g["d"] = a[dkey(a)] + f * (b[dkey(b)] - a[dkey(a)])
            return g
    return None


def load(case_dir="."):
    """The case's classes and parcels, every number the gates and the panels use."""
    p = lambda *a: os.path.join(case_dir, *a)  # noqa: E731
    ini = read_ini(p("input.ini"))
    lx, ys, zs, gas = read_solfile(p("INPUT", "solfile.tec"))
    tmin, tmax, cols = load_properties(p("INPUT", "properties.dat"))
    case = {"model": ini[("IGLOO-Models", "breakup")], "sigma": float(ini[("IGLOO-Properties", "sigma")]),
            "mul": float(ini[("IGLOO-Properties", "mu")]), "rhog": gas["Roi( 1 )"], "ug": gas["U"],
            "mug": gas["MIL"], "tg": gas["T"], "lx": lx}
    case["rhol"] = table_value(cols["rho"], tmin, tmax, case["tg"])
    classes = []
    for m, n, kv, rp in sorted(read_bc(p("INPUT", "bc.txt")), key=lambda r: (r[1], r[0])):
        c = {"m": m, "n": n, "kv": kv, "d0": 2.0 * rp, "u0": (1.0 - kv) * case["ug"],
             "y": (ys[m - 1], ys[m]), "z": (zs[n - 1], zs[n]), "ids": []}
        c["we"] = case["rhog"] * c["d0"] * c["u0"] ** 2 / case["sigma"]
        c["re"] = case["rhog"] * c["d0"] * c["u0"] / case["mug"]
        c["tstar"] = c["d0"] * math.sqrt(case["rhol"] / case["rhog"]) / c["u0"]
        c["tb"] = TB_STAR * c["tstar"]
        c["hf92"] = hf92_smd(c, case)
        classes.append(c)
    case["classes"] = classes
    case["oh"] = case["mul"] / math.sqrt(case["rhol"] * classes[0]["d0"] * case["sigma"])

    traj, exits = {}, {}
    for r in rows(p("OUTPUT", "trajectories-A.dat")):
        traj.setdefault(r["ID"], []).append(r)
    for r in rows(p("OUTPUT", "outloc-A.dat")):
        exits.setdefault(r["ID"], []).append(r)
    parcels, orphans = {}, []
    for pid in sorted(set(traj) | set(exits)):
        seq = sorted(traj.get(pid, []), key=lambda r: (r["t"], r["X"], r["Y"], r["Z"], r["U"], r["V"], r["W"],
                                                        r[dkey(r)]))
        ex = exits.get(pid, [])
        if not seq or len(ex) != 1 or seq[0]["t"] != 0.0:
            orphans.append((pid, len(seq), len(ex)))
            continue
        inj, last, ex = seq[0], seq[-1], ex[0]
        cand = [c for c in classes if c["y"][0] - X_TOL <= inj["Y"] <= c["y"][1] + X_TOL
                and c["z"][0] - X_TOL <= inj["Z"] <= c["z"][1] + X_TOL]
        hit = [c for c in cand if abs(speed(inj) - c["kv"] * case["ug"]) <= U_TOL]
        if len(hit) != 1:
            orphans.append((pid, len(seq), 1))
            continue
        c = hit[0]
        c["ids"].append(pid)
        dk = dkey(seq[0])
        half = [r[dk] for r in seq if r["X"] >= 0.5 * lx]
        du = (case["ug"] - last["U"], -last["V"], -last["W"])
        tb = at_time(seq, c["tb"])
        parcels[pid] = {
            "cls": classes.index(c), "d_inj": inj[dk], "d": ex[dkey(ex)], "mdot": ex["mdot"], "x_exit": ex["X"],
            "bc": ex["bc"], "events": sum(1 for a, b in zip(seq, seq[1:]) if b[dk] != a[dk]),
            "half_const": bool(half) and max(half) == min(half),
            "we_last": case["rhog"] * sum(v * v for v in du) * 0.5 * last[dk] / case["sigma"],
            "t_last": last["t"], "tb": tb,
        }
        if tb is not None:
            ub = math.sqrt((case["ug"] - tb["U"]) ** 2 + tb["V"] ** 2 + tb["W"] ** 2)
            parcels[pid]["ub"] = ub
    case["parcels"], case["orphans"] = parcels, orphans
    for c in classes:
        ps = [parcels[i] for i in c["ids"]]
        if ps:
            c["smd"] = sum(q["mdot"] for q in ps) / sum(q["mdot"] / q["d"] for q in ps)
            c["mmd"] = mmd([q["d"] for q in ps], [q["mdot"] for q in ps])
    return case


def table(case):
    print("%s products, d0 %.0f um, rho_L/rho_G %.1f, Oh %.5f, U %.0f m/s, box %.2f m"
          % (case["model"], case["classes"][0]["d0"] * 1e6, case["rhol"] / case["rhog"], case["oh"], case["ug"],
             case["lx"]))
    print(" cls     We    u0[m/s]    Re   n  SMD/d0 [HF92]14  MMD/SMD events lat   x(5t*)[m]   We_r,last")
    for k, c in enumerate(case["classes"]):
        ps = [case["parcels"][i] for i in c["ids"]]
        if not ps:
            print(" %3d %7.2f %8.3f %6.0f   0" % (k + 1, c["we"], c["u0"], c["re"]))
            continue
        xb = [q["tb"]["X"] for q in ps if q["tb"] is not None]
        ev = [q["events"] for q in ps]
        lat = sum(1 for q in ps if q["x_exit"] < case["lx"] - X_TOL)
        print(" %3d %7.2f %8.3f %6.0f %3d  %6.4f  %6.4f   %6.3f  %2d-%-2d %3d  %s  %6.3f"
              % (k + 1, c["we"], c["u0"], c["re"], len(ps), c["smd"] / c["d0"], c["hf92"], c["mmd"] / c["smd"],
                 min(ev), max(ev), lat, ("%5.3f-%5.3f" % (min(xb), max(xb))) if xb else "    -      ",
                 max(q["we_last"] for q in ps)))


def gates(case):
    fails = 0
    cl, ps, lx = case["classes"], case["parcels"], case["lx"]
    counts = [len(c["ids"]) for c in cl]
    ok = not case["orphans"] and min(counts) >= N_MIN
    print("[%s] P0 population: %d classes, %d parcels, %d..%d per class (>= %d), %d parcel(s) without a class or "
          "one exit record%s" % ("PASS" if ok else "FAIL", len(cl), len(ps), min(counts), max(counts), N_MIN,
                                 len(case["orphans"]),
                                 "" if not case["orphans"] else " (first: ID %d, %d rows, %d exits)"
                                 % case["orphans"][0]))
    fails += not ok
    r = max(q["d"] / cl[q["cls"]]["d0"] for q in ps.values()) if ps else float("nan")
    ok = bool(ps) and r < BROKE
    print("[%s] P1 every drop broke: max d_exit/d0 = %.6f (< %.3f)" % ("PASS" if ok else "FAIL", r, BROKE))
    fails += not ok
    lat = sum(1 for q in ps.values() if q["x_exit"] < lx - X_TOL)
    var = sum(1 for q in ps.values() if not q["half_const"])
    wer = max((q["we_last"] for q in ps.values()), default=float("nan"))
    ok = bool(ps) and lat == 0 and var == 0 and wer < WE_R_MAX
    print("[%s] P2 the cascade ended inside the box: %d of %d parcels leave before the outlet, %d change d over "
          "x >= %g m, max We_r on the last row %.3f (< %.0f)"
          % ("PASS" if ok else "FAIL", lat, len(ps), var, 0.5 * lx, wer, WE_R_MAX))
    fails += not ok
    gap = [q["t_last"] / cl[q["cls"]]["tb"] for q in ps.values()]
    nb = sum(1 for q in ps.values() if q["tb"] is None)
    ok = bool(ps) and nb == 0
    print("[%s] P3 t_b = 5 t* bracketed: %d of %d parcels leave before t_b, min t_last/t_b = %.3f"
          % ("PASS" if ok else "FAIL", nb, len(ps), min(gap) if gap else float("nan")))
    fails += not ok
    return fails


def main(case_dir="."):
    case = load(case_dir)
    table(case)
    fails = gates(case)
    print("%s: %s" % (os.path.basename(os.path.abspath(case_dir)), "PASS" if not fails else "%d gate(s) FAIL" % fails))
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main())
