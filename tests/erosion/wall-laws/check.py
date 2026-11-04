#!/usr/bin/env python3
"""
wall-laws -- the erosion laws of [IGLOO-Erosion] in OUTPUT/wall.tec: Oka et al. (2005) for material A,
Neilson-Gilchrist (1968) for material B, on wall-flux's box (face 2, x = 0.15, coded 301; A in two groups, B in one).

Eight assigned NoDrag streams hit face 2 on straight lines at their injection velocities. For each, this oracle
computes the impact angle from the surface, alpha = asin(u/|v|), and the two laws as the papers write them:
  Oka   E = g(alpha) E90 [mm^3/kg], g = sin(alpha)^n1 (1 + Hv (1 - sin(alpha)))^n2,
        E90 = K (a Hv)^(k1 b) (v/v')^k2 (D/D')^k3; er = 1e-9 rho_wall E [kg/kg]
        (Part 1 eqs. 1-2, Part 2 eq. 4; n1, n2, k2 = s Hv^q from the deck's s and q)
  N-G   E = V^2 cos^2(alpha) sin(n alpha)/(2 phi) below alpha0 = pi/(2n), V^2 cos^2(alpha)/(2 phi) above it,
        plus (V sin(alpha) - k)^2/(2 epsilon) where V sin(alpha) > k (eqs. 3-5; the deck's k is the paper's K)

Gates
  L0  the run ends normally; the log names A: Oka and B: Neilson-Gilchrist with the constants it uses, n1, n2 and
      k2 the correlations s Hv^q (ES12.5 prints)
  L1  er-A and er-B per cell = sum(mdot E)/sum(mdot) over the streams that hit it (the groups cancel), 0 elsewhere
  L2  g-A, g-B per cell = groups x sum(mdot)/(dy dz); rdot = (g-A er-A + g-B er-B)/rho_wall
  L3  the deck covers both branches of Neilson-Gilchrist and its threshold: a stream with alpha < alpha0 and
      V sin(alpha) <= k, one with alpha > alpha0 and V sin(alpha) > k; the two streams of cell (1,1) differ in angle
  L4  ref-direct/: this deck with n1, n2 and k2 given as values (the correlations, 17 digits): the same wall.tec
"""
import configparser
import math
import os
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "tools"))
import walltec  # noqa: E402

IGLOO = os.environ.get("IGLOO_EXE", os.path.abspath(os.path.join(HERE, "..", "..", "..", "bin", "IGLOO")))
LX, LY, LZ = 0.15, 0.05, 0.05
NY, NZ = 5, 5
DY, DZ = LY / NY, LZ / NZ
GROUPS = {"A": 2, "B": 1}
REL = 1.0e-12                 # 16-digit prints of sums of a few products
REL_LOG = 6.0e-6              # ES12.5 prints


def fail(msg):
    print(f"[FAIL] {msg}")
    return 1


def read_deck(path):
    cfg = configparser.ConfigParser(inline_comment_prefixes=(";",))
    cfg.optionxform = str
    cfg.read(path)
    return cfg


def floats(cfg, sec, key):
    return [float(v) for v in cfg.get(sec, key).split()]


def oka(v, alpha, d, hv, a, b, K, k1, k2, k3, n1, n2, vref, dref):
    s = math.sin(alpha)
    g = s ** n1 * (1.0 + hv * (1.0 - s)) ** n2
    e90 = K * (a * hv) ** (k1 * b) * (v / vref) ** k2 * (d / dref) ** k3
    return g * e90


def ng(v, alpha, phi, eps, k, n):
    a0 = math.pi / (2.0 * n)
    cut = v * v * math.cos(alpha) ** 2 * (math.sin(n * alpha) if alpha < a0 else 1.0) / (2.0 * phi)
    vn = v * math.sin(alpha) - k
    return cut + (vn * vn / (2.0 * eps) if vn > 0.0 else 0.0)


def constants(cfg):
    """Material A's Oka constants (n1, n2, k2 from the deck's values or s, q) and B's Neilson-Gilchrist ones."""
    E = "IGLOO-Erosion"
    hv, a, b = (float(cfg.get(E, k)) for k in ("Hv", "a", "b"))

    def form(key):
        if cfg.has_option(E, key):
            return floats(cfg, E, key)[0]
        return floats(cfg, E, key + "-s")[0] * hv ** floats(cfg, E, key + "-q")[0]

    A = dict(hv=hv, a=a, b=b, K=floats(cfg, E, "K")[0], k1=floats(cfg, E, "k1")[0], k2=form("k2"),
             k3=floats(cfg, E, "k3")[0], n1=form("n1"), n2=form("n2"), vref=floats(cfg, E, "v-ref")[0],
             dref=floats(cfg, E, "D-ref")[0])
    B = dict(phi=floats(cfg, E, "phi")[1], eps=floats(cfg, E, "epsilon")[1], k=floats(cfg, E, "k")[1],
             n=floats(cfg, E, "n")[1])
    return A, B, float(cfg.get(E, "wall-density"))


def streams(cfg):
    """(cell, alpha, |v|, d, mdot) of every stream, on its straight line to x = LX."""
    S = "IGLOO-BC"
    x, y, z = floats(cfg, S, "x"), floats(cfg, S, "y"), floats(cfg, S, "z")
    up, vp, wp = floats(cfg, S, "up"), floats(cfg, S, "vp"), floats(cfg, S, "wp")
    diam, mdot = floats(cfg, S, "diam"), floats(cfg, S, "mdot")
    out = []
    for s in range(len(x)):
        t = (LX - x[s]) / up[s]
        yw, zw = y[s] + vp[s] * t, z[s] + wp[s] * t
        speed = math.sqrt(up[s] ** 2 + vp[s] ** 2 + wp[s] ** 2)
        out.append(((int(yw // DY) + 1, int(zw // DZ) + 1), math.asin(up[s] / speed), speed, diam[s], mdot[s]))
    return out


def expected(cfg):
    A, B, rho = constants(cfg)
    sums = {}
    for cell, alpha, v, d, m in streams(cfg):
        eA = 1.0e-9 * rho * oka(v, alpha, d, A["hv"], A["a"], A["b"], A["K"], A["k1"], A["k2"], A["k3"], A["n1"],
                                A["n2"], A["vref"], A["dref"])
        eB = ng(v, alpha, B["phi"], B["eps"], B["k"], B["n"])
        acc = sums.setdefault(cell, [0.0, 0.0, 0.0])
        acc[0] += m * eA
        acc[1] += m * eB
        acc[2] += m
    exp = {}
    for cell, (mA, mB, m) in sums.items():
        gA, gB = GROUPS["A"] * m / (DY * DZ), GROUPS["B"] * m / (DY * DZ)
        exp[cell] = {"er-A": mA / m, "er-B": mB / m, "g-A": gA, "g-B": gB,
                     "rdot": (gA * mA / m + gB * mB / m) / rho}
    return exp, A, B


def cell_values(zone, name):
    vals = zone["cells"][name]
    return {(j + 1, k + 1): vals[j + NY * k] for k in range(NZ) for j in range(NY)}


def compare(zone, exp, label):
    """Worst relative difference of every variable against the oracle; exact zero where no stream hits."""
    rc, worst = 0, {}
    for name in ("er-A", "er-B", "g-A", "g-B", "rdot"):
        got = cell_values(zone, name)
        w = 0.0
        for c, val in got.items():
            if c in exp:
                w = max(w, abs(val - exp[c][name]) / abs(exp[c][name]))
            elif val != 0.0:
                rc |= fail(f"{label} {name}: {val!r} in cell {c}, which no stream hits")
        worst[name] = w
    return rc, worst


def run_direct(cfg, A):
    """This deck with n1, n2 and k2 as values, into ref-direct/."""
    ref = os.path.join(HERE, "ref-direct")
    shutil.rmtree(ref, ignore_errors=True)
    os.makedirs(os.path.join(ref, "OUTPUT"))
    os.symlink(os.path.join(HERE, "INPUT"), os.path.join(ref, "INPUT"))
    lines = []
    with open(os.path.join(HERE, "input.ini")) as fh:
        for ln in fh:
            key = ln.split("=")[0].strip()
            if key in ("k2-s", "k2-q", "n1-s", "n1-q", "n2-s", "n2-q"):
                continue
            lines.append(ln)
            if key == "k1":
                for k in ("k2", "n1", "n2"):
                    lines.append(f"{k:<12} = {A[k]!r} 0.0\n")
    with open(os.path.join(ref, "input.ini"), "w") as fh:
        fh.writelines(lines)
    with open(os.path.join(ref, "run_out.txt"), "w") as out, open(os.path.join(ref, "run_err.txt"), "w") as err:
        r = subprocess.run([IGLOO], cwd=ref, stdout=out, stderr=err)
    return ref, r.returncode


def logged(log, prefix):
    """The constants line after the log line `prefix`: {name: value}."""
    lines = log.splitlines()
    for i, ln in enumerate(lines):
        if ln == prefix and i + 1 < len(lines):
            out = {}
            for item in lines[i + 1].split(","):
                k, v = item.split()
                out[k] = float(v)
            return out
    return None


def main():
    rc = 0
    cfg = read_deck(os.path.join(HERE, "input.ini"))
    exp, A, B = expected(cfg)
    log = open(os.path.join(HERE, "run_out.txt")).read()
    err = open(os.path.join(HERE, "run_err.txt")).read()

    # L0 -------------------------------------------------------------------------------------------------------------
    la = logged(log, "         material A: Oka")
    lb = logged(log, "         material B: Neilson-Gilchrist")
    wantA = {"K": A["K"], "k1": A["k1"], "k2": A["k2"], "k3": A["k3"], "n1": A["n1"], "n2": A["n2"],
             "v-ref": A["vref"], "D-ref": A["dref"]}
    wantB = {"phi": B["phi"], "epsilon": B["eps"], "k": B["k"], "n": B["n"]}
    bad = []
    if "Stop condition" not in log or err.strip() or la is None or lb is None:
        bad.append("the run or its model lines")
    else:
        for want, got in ((wantA, la), (wantB, lb)):
            for k, v in want.items():
                if k not in got or abs(got[k] - v) > REL_LOG * abs(v):
                    bad.append(f"{k} {got.get(k)} for {v}")
    if bad:
        rc |= fail("L0 " + "; ".join(bad))
    else:
        print(f"L0 the log names A: Oka (n1 {la['n1']:.5e}, n2 {la['n2']:.5e}, k2 {la['k2']:.5e} = s Hv^q) and "
              "B: Neilson-Gilchrist with the deck's constants PASS")

    # L1, L2 ---------------------------------------------------------------------------------------------------------
    wall = walltec.read(os.path.join(HERE, "OUTPUT", "wall.tec"))
    zone = wall["zones"][0]
    r, worst = compare(zone, exp, "L1")
    rc |= r
    for name in ("er-A", "er-B"):
        if worst[name] > REL:
            rc |= fail(f"L1 {name}: worst {worst[name]:.2e} relative")
        else:
            print(f"L1 {name} in the {len(exp)} hit cells = sum(mdot E)/sum(mdot) of the oracle: worst "
                  f"{worst[name]:.1e} <= {REL:.0e}, 0 elsewhere PASS")
    for name in ("g-A", "g-B", "rdot"):
        if worst[name] > REL:
            rc |= fail(f"L2 {name}: worst {worst[name]:.2e} relative")
        else:
            print(f"L2 {name}: worst {worst[name]:.1e} <= {REL:.0e} PASS")

    # L3 -------------------------------------------------------------------------------------------------------------
    st = streams(cfg)
    a0 = math.pi / (2.0 * B["n"])
    below = [i + 1 for i, s in enumerate(st) if s[1] < a0 and s[2] * math.sin(s[1]) <= B["k"]]
    above = [i + 1 for i, s in enumerate(st) if s[1] > a0 and s[2] * math.sin(s[1]) > B["k"]]
    shared = [s[1] for s in st if s[0] == (1, 1)]
    if not below or not above or len(shared) != 2 or abs(shared[0] - shared[1]) < 1e-3:
        rc |= fail(f"L3 coverage: below alpha0 and k {below}, above {above}, cell (1,1) angles {shared}")
    else:
        print(f"L3 streams {below} below alpha0 = {math.degrees(a0):.0f} deg and k, {above} above; cell (1,1) at "
              f"{math.degrees(shared[0]):.0f} and {math.degrees(shared[1]):.0f} deg PASS")

    # L4 -------------------------------------------------------------------------------------------------------------
    ref, code = run_direct(cfg, A)
    if code != 0:
        return rc | fail(f"L4 the run with n1, n2, k2 as values ended with {code}")
    zd = walltec.read(os.path.join(ref, "OUTPUT", "wall.tec"))["zones"][0]
    w = 0.0
    for name in ("er-A", "er-B", "rdot"):
        a, b = cell_values(zone, name), cell_values(zd, name)
        for c in a:
            if a[c] != 0.0 or b[c] != 0.0:
                w = max(w, abs(a[c] - b[c]) / max(abs(a[c]), abs(b[c])))
    if w > 1.0e-13:
        rc |= fail(f"L4 ref-direct/ differs from the correlations' run: worst {w:.2e} relative")
    else:
        print(f"L4 n1, n2, k2 as values (ref-direct/) give the same er-A, er-B, rdot: worst {w:.1e} <= 1e-13 PASS")

    if rc:
        print("[FAIL] wall-laws violation(s) -- see above.")
    return rc


if __name__ == "__main__":
    sys.exit(main())
