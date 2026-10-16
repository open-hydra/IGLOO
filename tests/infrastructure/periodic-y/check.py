#!/usr/bin/env python3
"""Independent oracle for the TRANSLATIONAL PERIODIC (bcdef 201) e2e case.

Physics is the standard/body-force decoupling (kV=kT=1, Stokes exact-linear drag,
transverse g_y) with g_y = -8000 so each particle drifts ~2.8 box heights while
crossing the domain; faces 3/4 (y-min/y-max) are a 201 periodic pair.

What 201 must do (obj_particles.updateCell -> periodicTransport): translate the
exiting particle by T = partner_face_center - exit_face_center = +-Ly*yhat with
VELOCITY UNCHANGED, then relocate it in the block. Gates:

  1. v_y(x) follows the anchored closed form v = v_inf + (v_a - v_inf)e^{-(x-x_a)/L}
     across ALL wraps -- velocity mangling at any transport breaks the exponential;
  2. y(x) equals the unwrapped closed form MODULO Ly (the transport is an exact
     +-Ly translation; any drift or reflection shows up as a modular residual);
  3. every particle wraps >= MIN_WRAPS times (non-vacuous) and exits at x = Lx;
  4. u == u_g, w == 0, T == T_g decoupling guards hold through the transports;
  5. (S-P) the field mollifier crosses the periodic faces: a run with `mollify = off` smoothed here
     (8 passes of the [1/4, 1/2, 1/4] filter, zero flux in x and z, periodic in y, the box's uniform
     volumes) equals the smoothed run's source.tec field by field to 1e-12 of the field's scale --
     (a) this case against its unsmoothed sibling (ref-raw/); (b) the pair again with only face 1's
     first row injecting (ref-row/, ref-row-raw/). The case's parcels start one cell apart, so its
     deposits are uniform in y to ~1e-6; (b)'s are not, and there the same passes with zero flux in y
     must differ by more than 1e-6 (the periodic faces carry a flux the gate can see).
"""
import math
import os
import re
import shutil
import subprocess
import sys

TRAJ   = "OUTPUT/trajectories-A.dat"
OUTLOC = "OUTPUT/outloc-A.dat"
SOURCE = "OUTPUT/source.tec"
IGLOO  = os.environ.get("IGLOO_EXE",
                        os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", "bin", "IGLOO")))

# ---- known test inputs (SI). NOT read from production output. ----------------
RHO_P = 2950.0
D_P   = 1.189e-5
MU    = 1.8e-5
U_G   = 10.0
G_Y   = -8000.0
LY    = 0.05
LX    = 0.15
T_G   = 600.0

TAU   = RHO_P * D_P**2 / (18.0 * MU)     # 1.287e-3 s
L     = U_G * TAU                        # 12.87 mm
V_INF = G_Y * TAU                        # -10.30 m/s

# ---- error model --------------------------------------------------------------
HALF_ULP  = 0.5e-6
INT_FLOOR = 1.0e-8
V_BAND    = 1.0e-3      # exclude the v ~= v_inf saturation band (F12.6-unresolved)
TOL_Y     = 5.0e-6      # modular y residual [m]: measured 5e-7 (anchor truncation + print ULPs)
MIN_PTS   = 5
N_GOOD    = 20
MIN_WRAPS = 2
TOL_U     = 1.0e-4
TOL_W     = 1.0e-6
TOL_T_REL = 1.0e-3
X_EXIT    = 0.14
THETA     = 0.25        # the mollifier's per-sweep coefficient
NPASS     = 8           # its default number of passes
TOL_FIELD = 1.0e-12     # scale-relative, compare_tec.py's rule
MIN_SEEN  = 1.0e-6      # periodic vs zero-flux smoothing: the gate's non-vacuity floor
MOLLIFY_LINE = "Field mollification ON: 8 binomial pass(es)"


def load(path):
    parts = {}
    with open(path) as f:
        for line in f:
            c = line.split()
            if len(c) < 10:
                continue
            try:
                x, y, z, u, v, w, T, dp = (float(c[i]) for i in range(8))
                pid = int(c[9])
            except ValueError:
                continue
            parts.setdefault(pid, []).append((x, y, u, v, w, T))
    return parts


def read_cells(path):
    """One-zone BLOCK-packed source.tec: cells (ni, nj, nk) and {cell variable: values, i fastest}."""
    txt = open(path).read()
    names = re.search(r"VARIABLES\s*=(.*?)\n", txt).group(1).split()
    head, body = re.split(r"\n\s*ZONE", txt, maxsplit=1)[1].split("\n", 1)
    I, J, K = (int(re.search(r"\b%s\s*=\s*(\d+)" % c, head).group(1)) for c in "IJK")
    a, b = (int(x) for x in re.search(r"\[(\d+)-(\d+)\]\s*=\s*CELLCENTERED", head).groups())
    vals = [float(v) for v in body.split()]
    nnode, ncell = I * J * K, (I - 1) * (J - 1) * (K - 1)
    cells, pos = {}, 0
    for k, name in enumerate(names, 1):
        n = ncell if a <= k <= b else nnode
        if a <= k <= b:
            cells[name] = vals[pos:pos + n]
        pos += n
    if pos != len(vals):
        raise ValueError(f"{path}: {len(vals)} values, the layout needs {pos}")
    return I - 1, J - 1, K - 1, cells


def smooth(q, ni, nj, nk, periodic):
    """NPASS passes of the conservative binomial filter on uniform volumes, sweeps in x, y and (3D) z;
    zero flux at every face but y's when periodic (the first and last rows are then neighbours)."""
    q = list(q)
    idx = lambda i, j, k: (k * nj + j) * ni + i
    dirs = [(1, 0, 0), (0, 1, 0)] + ([(0, 0, 1)] if nk > 1 else [])
    for _ in range(NPASS):
        for d in dirs:
            new = [0.0] * len(q)
            for k in range(nk):
                for j in range(nj):
                    for i in range(ni):
                        c, phi = idx(i, j, k), 0.0
                        for s in (1, -1):
                            ii, jj, kk = i + s * d[0], j + s * d[1], k + s * d[2]
                            if d[1] and periodic:
                                jj %= nj
                            if 0 <= ii < ni and 0 <= jj < nj and 0 <= kk < nk:
                                phi += THETA * (q[idx(ii, jj, kk)] - q[c])
                        new[c] = q[c] + phi
            q = new
    return q


def row_bc(src, dst):
    """The case's bc.txt with every face-1 inlet record but the first row's (j = 1) rewritten as a
    header-only 400; returns the number of 401 records kept. Box-fixture headers: b face m n 1 code."""
    lines = open(src).read().split("\n")
    out, i, kept = [], 0, 0
    while i < len(lines):
        t = lines[i].split()
        if len(t) != 6:
            out.append(lines[i]); i += 1
            continue
        code = int(t[5])
        npay = 1 if code in (101, 201, 401, 402, 403) else 0
        if code == 401 and int(t[1]) == 1 and int(t[2]) != 1:
            out.append("   " + "".join(f"{int(x):6d}" for x in t[:5]) + f"{400:6d}")
        else:
            kept += code == 401
            out.extend(lines[i:i + 1 + npay])
        i += 1 + npay
    open(dst, "w").write("\n".join(out))
    return kept


def sibling(name, extra_ini, row_only):
    """This case run again into <name>/ with extra_ini in [IGLOO-General]; with row_only only face 1's
    first row injects. Returns (dir, exit code)."""
    case = os.getcwd()
    d = os.path.join(case, name)
    shutil.rmtree(d, ignore_errors=True)
    os.makedirs(os.path.join(d, "OUTPUT"))
    if row_only:
        os.makedirs(os.path.join(d, "INPUT"))
        for f in os.listdir(os.path.join(case, "INPUT")):
            if f != "bc.txt":
                os.symlink(os.path.join(case, "INPUT", f), os.path.join(d, "INPUT", f))
        row_bc(os.path.join(case, "INPUT", "bc.txt"), os.path.join(d, "INPUT", "bc.txt"))
    else:
        os.symlink(os.path.join(case, "INPUT"), os.path.join(d, "INPUT"))
    ini = open(os.path.join(case, "input.ini")).read()
    if extra_ini:
        ini = ini.replace("[IGLOO-General]\n", "[IGLOO-General]\n" + extra_ini + "\n", 1)
    open(os.path.join(d, "input.ini"), "w").write(ini)
    with open(os.path.join(d, "run_out.txt"), "w") as out, open(os.path.join(d, "run_err.txt"), "w") as err:
        code = subprocess.run([IGLOO], cwd=d, stdout=out, stderr=err).returncode
    return d, code


def compare_smoothed(lab, smoothed_dir, raw_dir, need_seen):
    """The smoothed run's source.tec against the raw run's smoothed here with the periodic pair joined,
    field by field; with need_seen the zero-flux smoothing must differ by more than MIN_SEEN."""
    if MOLLIFY_LINE not in open(os.path.join(smoothed_dir, "run_out.txt")).read():
        print(f"[FAIL] S-P {lab}: run_out.txt does not report '{MOLLIFY_LINE}'")
        return 1
    ni, nj, nk, got = read_cells(os.path.join(smoothed_dir, SOURCE))
    rni, rnj, rnk, raw = read_cells(os.path.join(raw_dir, SOURCE))
    if (ni, nj, nk) != (rni, rnj, rnk) or set(got) != set(raw):
        print(f"[FAIL] S-P {lab}: the unsmoothed run's source.tec has another layout")
        return 1
    rc, seen = 0, 0.0
    for name in got:
        per = smooth(raw[name], ni, nj, nk, True)
        zf = smooth(raw[name], ni, nj, nk, False)
        scale = max(abs(v) for v in got[name])
        if scale == 0.0:
            ok = all(v == 0.0 for v in per)
            print(f"S-P {lab} {name}: zero field, oracle zero  {'PASS' if ok else 'FAIL'}")
            rc |= 0 if ok else 1
            continue
        err = max(abs(a - b) for a, b in zip(got[name], per)) / scale
        dzf = max(abs(a - b) for a, b in zip(zf, per)) / scale
        seen = max(seen, dzf)
        verdict = "PASS" if err <= TOL_FIELD else "FAIL"
        print(f"S-P {lab} {name}: against the periodic smoothing {err:.3e} of the field scale (tol {TOL_FIELD:.0e}); "
              f"zero-flux smoothing differs by {dzf:.3e}  {verdict}")
        rc |= 0 if err <= TOL_FIELD else 1
    if need_seen and seen <= MIN_SEEN:
        print(f"[FAIL] S-P {lab}: periodic and zero-flux smoothing differ by {seen:.3e} <= {MIN_SEEN:.0e} -- blind")
        rc = 1
    return rc


def smoothing_gate():
    """S-P: (a) this run against its unsmoothed sibling (ref-raw/) smoothed here; (b) the same pair with
    only face 1's first row injecting (ref-row/, ref-row-raw/), whose deposits are not uniform in y. The
    case's own deposits are uniform in y to ~1e-6 (its parcels start one cell apart): (a) checks the
    artefact itself, (b) carries the non-vacuity floor."""
    rc = 0
    raw, code = sibling("ref-raw", "mollify = off", False)
    if code != 0:
        print(f"[FAIL] S-P the unsmoothed run exited {code}")
        return 1
    rc |= compare_smoothed("case", ".", raw, False)
    row, code = sibling("ref-row", "", True)
    rowraw, code2 = sibling("ref-row-raw", "mollify = off", True)
    if code != 0 or code2 != 0:
        print(f"[FAIL] S-P the first-row runs exited {code}, {code2}")
        return 1
    rc |= compare_smoothed("first row", row, rowraw, True)
    return rc


def main():
    try:
        parts = load(TRAJ)
    except FileNotFoundError:
        print(f"[FAIL] {TRAJ} not found -- did the solver run?")
        return 1
    if not parts:
        print(f"[FAIL] no trajectory rows in {TRAJ}")
        return 1
    print(f"oracle: tau={TAU:.6e} s  L={L:.4e} m  v_inf={V_INF:.4f} m/s  "
          f"drift/height={abs(V_INF)*(LX/U_G - TAU)/LY:.2f} wraps expected\n")

    nviol = 0
    verified = 0
    worst_v = (0.0, "")
    worst_y = (0.0, "")
    total_wraps = {}

    for pid in sorted(parts):
        rows = sorted(parts[pid], key=lambda r: r[0])
        if len(rows) <= 1:
            continue
        xa, ya, _, va, _, _ = rows[0]
        npts = 0
        wraps = 0
        y_prev = ya
        bad = False
        for (x, y, u, v, w, T) in rows:
            if abs(u - U_G) > TOL_U:
                print(f"[FAIL] ID={pid}: u={u:.6f} != u_g -- decoupling broken")
                nviol += 1; bad = True; break
            if abs(w) > TOL_W:
                print(f"[FAIL] ID={pid}: w={w:.2e} != 0")
                nviol += 1; bad = True; break
            if abs(T - T_G) > TOL_T_REL * T_G:
                print(f"[FAIL] ID={pid}: T={T:.3f} != T_g")
                nviol += 1; bad = True; break

            dt = (x - xa) / U_G
            # wrap counter: consecutive-row y jump ~ Ly that is NOT physical motion
            if abs(y - y_prev) > 0.5 * LY:
                wraps += 1
            y_prev = y

            # gate 1: velocity closed form (transport must not touch v)
            v_pred = V_INF + (va - V_INF) * math.exp(-dt / TAU)
            if abs(v - V_INF) > V_BAND:
                tol = HALF_ULP * (2.0 + 2.0 * abs(v - V_INF) / L) + INT_FLOOR
                resid = abs(v - v_pred)
                npts += 1
                if resid / tol > worst_v[0]:
                    worst_v = (resid / tol, f"ID={pid} x={x:.4f}")
                if resid > tol:
                    print(f"[FAIL] ID={pid} x={x:.4f}: |v-v_pred|={resid:.3e} > {tol:.3e}")
                    nviol += 1; bad = True

            # gate 2: position closed form modulo Ly (transport = exact +-Ly)
            y_unwrapped = ya + V_INF * dt + (va - V_INF) * TAU * (1.0 - math.exp(-dt / TAU))
            r = (y - y_unwrapped) % LY
            resid_y = min(r, LY - r)
            if resid_y / TOL_Y > worst_y[0]:
                worst_y = (resid_y / TOL_Y, f"ID={pid} x={x:.4f}")
            if resid_y > TOL_Y:
                print(f"[FAIL] ID={pid} x={x:.4f}: modular y residual "
                      f"{resid_y:.3e} > {TOL_Y:.1e}")
                nviol += 1; bad = True
        total_wraps[pid] = wraps
        if not bad and npts >= MIN_PTS and wraps >= MIN_WRAPS:
            verified += 1

    exits = 0
    try:
        for ln in open(OUTLOC).read().splitlines()[2:]:
            c = ln.split()
            if len(c) == 9 and float(c[0]) > X_EXIT:
                exits += 1
    except FileNotFoundError:
        print(f"[FAIL] {OUTLOC} not found")
        return 1

    wr = sorted(total_wraps.values())
    print(f"particles: {len(parts)}  outlet exits (x>{X_EXIT}): {exits}")
    print(f"wraps per particle: min={wr[0]} max={wr[-1]}  "
          f"(need >= {MIN_WRAPS} each for the gate)")
    print(f"worst v point: resid/tol={worst_v[0]:.3f} ({worst_v[1]})")
    print(f"worst modular-y point: resid/tol={worst_y[0]:.3f} ({worst_y[1]})")
    print(f"verified wrappers: {verified} (need >= {N_GOOD})")

    sp = smoothing_gate()
    ok = (nviol == 0 and verified >= N_GOOD and exits == len(parts) and sp == 0)
    print(("\n[PASS] 201 periodic transport verified: velocity untouched, "
           "position exact modulo Ly, all particles exit; the source smoothed across the periodic faces.") if ok
          else "\n[FAIL] periodic-transport or smoothing violation(s) -- see above.")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
