#!/usr/bin/env python3
"""temp-relax-varcp-euler -- the equivalent-Eulerian temperature of a material whose cp varies (ODE model 1).

temp-relax-varcp's run (enthalpy state on a varying-cp table spanning 250..800 K, four inlet groups at 210, 300,
750 and 840 K, kV = 1 so every parcel coasts at u_g along its own (j, k) column of cells) with out-file = ALL.
Model 1 deposits per cell the moments int |v| dt and int h |v| dt (h the enthalpy state); normalizeEUL divides
them and inverts the table, so with |v| constant the Eulerian temperature of cell c is

  T_p(c) = T(hbar_c),  hbar_c = (1/(t_out - t_in)) int_{t_in}^{t_out} h(T(t)) dt,

the time mean of the enthalpy over the parcel's residence in the cell. T(t) is temp-relax-varcp's closed form
(the table's segments marched from the injection state); on a segment of slope s, h(T(t)) is an exponential in t,
so the integral is exact, segment by segment. The residence window is read from the trajectory rows (out-time =
on): the row that enters the cell and the row that enters the next one (the exit record for the last cell).
The inverse T(h) is this check's own (the table's piecewise-linear enthalpy, the end segments extended).

Gates
  W   the witnesses: the enthalpy state, mollification off, the Eulerian field written
  E1  |T_p(c) - T(hbar_c)| <= TOL_E1 on every crossed cell
  E2  T_p(c) inside [T_in, T_out] of the two rows bounding the cell (within the F12.6 print)
  E3  coverage: all 1500 cells of the box (the 25 inlet columns) have rho_p > 0 and are checked
"""
import math
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "tools"))
sys.path.insert(0, os.path.join(HERE, "..", "temp-relax-varcp"))
import check as tr  # noqa: E402  temp-relax-varcp's inputs, table recipe and closed form
import datcols  # noqa: E402

EUL, TRAJ, EXIT, LOG = "OUTPUT/euler1.tec", "OUTPUT/trajectories-A.dat", "OUTPUT/outloc-A.dat", "run_out.txt"
NX, NY, NZ = 60, 5, 5
DX, DY, DZ = 0.15 / NX, tr.LY / NY, tr.LZ / NZ

TOL_E1 = 1.0e-5          # K: the moments are ODE states (rtol 1e-11), the window times print 9 digits
TOL_E2 = 1.0e-6          # K: the F12.6 print of the bounding rows
N_CELLS = NX * NY * NZ


def segment(T):
    """(i, s, T_b) of the table segment the relaxation from T is on, and the end T_b it runs to."""
    if T < tr.T_G:
        if T < tr.TMIN:
            return tr.TMIN, tr.slope(tr.TMIN), float(tr.TMIN)
        i = math.floor(T)
        return i, tr.slope(i), float(min(i + 1, tr.T_G))
    if T > tr.TMAX:
        return tr.TMAX - 1, tr.slope(tr.TMAX - 1), float(tr.TMAX)
    i = math.ceil(T) - 1
    return i, tr.slope(i), float(max(i, tr.T_G))


def h_integral(T0, t1, t2):
    """int_{t1}^{t2} h(T(t)) dt along the closed-form path from T(0) = T0, exact per table segment."""
    T, t, acc = T0, 0.0, 0.0
    for _ in range(4 * (tr.TMAX - tr.TMIN)):
        i, s, Tb = segment(T)
        rate = tr.K / s
        if Tb == tr.T_G or (T - tr.T_G) * (Tb - tr.T_G) <= 0.0:
            tb = math.inf
        else:
            tb = t + math.log((tr.T_G - T) / (tr.T_G - Tb)) / rate
        lo, hi = max(t, t1), min(tb, t2)
        if hi > lo:
            # h = H_i + s (T - i) = c0 + s (T - T_g), with T - T_g = (T_seg - T_g) exp(-rate (t - t_seg))
            c0 = tr.H[i] + s * (tr.T_G - i)
            acc += c0 * (hi - lo) + s * (T - tr.T_G) / rate * (math.exp(-rate * (lo - t)) - math.exp(-rate * (hi - t)))
        if tb >= t2:
            return acc
        T, t = Tb, tb
    raise RuntimeError("segment march did not converge")


def T_of_h(h):
    """Inverse of the table's piecewise-linear enthalpy, the end segments extended."""
    if h <= tr.H[tr.TMIN]:
        return tr.TMIN + (h - tr.H[tr.TMIN]) / tr.slope(tr.TMIN)
    if h >= tr.H[tr.TMAX - 1]:
        return (tr.TMAX - 1) + (h - tr.H[tr.TMAX - 1]) / tr.slope(tr.TMAX - 1)
    lo, hi = tr.TMIN, tr.TMAX - 1
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if tr.H[mid] <= h:
            lo = mid
        else:
            hi = mid
    return lo + (h - tr.H[lo]) / tr.slope(lo)


def read_euler(path):
    """BLOCK-packed Tecplot: the nodal x, y, z, then the cell-centred variables, i fastest."""
    lines = open(path).read().splitlines()
    names = re.findall(r"\S+", lines[0].split("=", 1)[1])
    zi = next(n for n, l in enumerate(lines) if l.strip().upper().startswith("ZONE"))
    I, J, K = (int(re.search(r"\b%s=(\d+)" % a, lines[zi]).group(1)) for a in "IJK")
    nn, nc = I * J * K, (I - 1) * (J - 1) * (K - 1)
    vals = [float(v) for v in " ".join(lines[zi + 1:]).split()]
    cells = {names[3 + a]: vals[3 * nn + a * nc: 3 * nn + (a + 1) * nc] for a in range(len(names) - 3)}
    return (I - 1, J - 1, K - 1), cells


def main():
    fails = []
    try:
        log = open(LOG, errors="replace").read()
        dims, cells = read_euler(EUL)
        traj = datcols.rows(TRAJ)
        exits = {r["ID"]: r for r in datcols.rows(EXIT)}
    except (OSError, ValueError, StopIteration, AttributeError) as exc:
        print("[FAIL] cannot read the outputs: %s" % exc)
        return 1
    for witness in ("Solving enthalpy equation", "Field mollification OFF"):
        if witness not in log:
            fails.append("W the log does not report '%s'" % witness)
    if dims != (NX, NY, NZ):
        fails.append("W the Eulerian field is %s cells, the box %s" % (dims, (NX, NY, NZ)))
    rho, tp = cells.get("rho<sub>p"), cells.get("T<sub>p")
    if rho is None or tp is None:
        print("[FAIL] euler1.tec carries no rho_p / T_p: %s" % sorted(cells))
        return 1

    parts = {}
    for r in traj:
        parts.setdefault(r["ID"], []).append(r)
    checked, worst1, worst2, seen = 0, (0.0, ""), (0.0, ""), set()
    for pid, rows in sorted(parts.items()):
        rows.sort(key=lambda r: r["t"])
        first = rows[0]
        j, k = int(first["Y"] / DY) + 1, int(first["Z"] / DZ) + 1
        T0 = tr.T_0[k]
        if first["t"] != 0.0 or abs(first["T"] - T0) > 1.0e-6:
            fails.append("W ID=%d: first row (t %.8e, T %.6f) is not the injection state (0, %.1f)"
                         % (pid, first["t"], first["T"], T0))
            continue
        ex = exits.get(pid)
        if ex is None:
            fails.append("E3 ID=%d: no exit record" % pid)
            continue
        bounds = [(r["X"], r["t"], r["T"]) for r in rows]
        if rows[-1]["t"] < ex["t"]:          # the trajectory ends at the exit: its last row is the exit record
            bounds.append((ex["X"], ex["t"], ex["T"]))
        for n in range(len(bounds) - 1):
            (xa, ta, Ta), (xb, tb, Tb) = bounds[n], bounds[n + 1]
            i = int(round(xa / DX)) + 1
            c = (i - 1) + NX * ((j - 1) + NY * (k - 1))
            seen.add(c)
            if not rho[c] > 0.0:
                fails.append("E3 cell (%d, %d, %d): rho_p = %g on a crossed cell" % (i, j, k, rho[c]))
                continue
            hbar = h_integral(T0, ta, tb) / (tb - ta)
            ref = T_of_h(hbar)
            r1 = abs(tp[c] - ref)
            if r1 / TOL_E1 > worst1[0]:
                worst1 = (r1 / TOL_E1, "cell (%d, %d, %d) T_p %.9f T(hbar) %.9f resid %.3e K"
                          % (i, j, k, tp[c], ref, r1))
            if r1 > TOL_E1:
                fails.append("E1 cell (%d, %d, %d): T_p %.9f vs T(hbar) %.9f (%.3e K)" % (i, j, k, tp[c], ref, r1))
            lo, hi = min(Ta, Tb) - TOL_E2, max(Ta, Tb) + TOL_E2
            out = max(lo - tp[c], tp[c] - hi, 0.0)
            if out / TOL_E2 > worst2[0]:
                worst2 = (out / TOL_E2, "cell (%d, %d, %d) T_p %.9f band [%.6f, %.6f]" % (i, j, k, tp[c], Ta, Tb))
            if out > 0.0:
                fails.append("E2 cell (%d, %d, %d): T_p %.9f outside [%.6f, %.6f]" % (i, j, k, tp[c], Ta, Tb))
            checked += 1
    others = [c for c in range(len(rho)) if rho[c] > 0.0 and c not in seen]
    if others:
        fails.append("E3 %d cell(s) with rho_p > 0 that no parcel crossed" % len(others))
    if checked != N_CELLS:
        fails.append("E3 %d cells checked, the 25 inlet columns hold %d" % (checked, N_CELLS))
    print("E1 worst: %s (%.3f of %.0e K)" % (worst1[1], worst1[0], TOL_E1))
    print("E2 worst: %s (%.3f of %.0e K outside)" % (worst2[1], worst2[0], TOL_E2))
    print("E3 %d cells checked of %d, %d parcels" % (checked, N_CELLS, len(parts)))
    for f in fails[:20]:
        print("[FAIL] " + f)
    if fails:
        print("\n[FAIL] %d violation(s)" % len(fails))
        return 1
    print("\n[PASS] the Eulerian temperature of every crossed cell is the table's T of the time-mean enthalpy")
    return 0


if __name__ == "__main__":
    sys.exit(main())
