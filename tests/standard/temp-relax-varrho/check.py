#!/usr/bin/env python3
"""Oracle for a parcel whose DENSITY varies with temperature: its diameter follows the density.

temp-relax-varcp's case (box, gas, the four inlet groups at 210/300/750/840 K, the varying cp; its check.py is
imported for the table recipe and the readers) with rho = 2950 - 0.5 (T - 250) kg/m^3 on 250..800 K. The
parcels exchange no mass, so each keeps the mass it was injected with, m = rho(T0) pi d0^3/6, and its diameter is

    d(T) = (6 m / (pi rho(T)))^(1/3)

Physics
-------
kV = 1 (Re = 0, Nu = 2, x = u_g t). With h(T) the table's piecewise-linear enthalpy (slope s on each segment,
the end segments extended):
    dh/dt = Nu k_g pi d(T) (T_g - T) / m,   so   t(T) = integral of s m / (Nu k_g pi d(T) (T_g - T)) dT
The reference T(x) inverts t(T) = x/u_g, integrating segment by segment with 20-point Gauss-Legendre (the
integrand is smooth on every segment inside the window below).

Gates
-----
V1  every row with |T - T_g| > T_BAND: the printed T against the reference, within the F12.6 truncation of T and
    of x (x/u_g is the time) plus the integrator floor
D1  every such row: the printed d against d(T_ref(x)), within its E13.6 half-ULP plus |dd/dT| times V1's budget
D2  every row: the printed d against d(T) at the printed T (half-ULP of both)
M   every row: the printed m equal to the group's injection mass (E13.6)
G   guards: u = u_g, v = w = 0, the first row at x = 0 and T0, the enthalpy-state and mollify-off witnesses,
    every parcel of every group verified (>= MIN_PTS window rows)
"""
import importlib.util
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "tools"))
from proptab import load_properties  # noqa: E402

_spec = importlib.util.spec_from_file_location("trv", os.path.join(HERE, "..", "temp-relax-varcp", "check.py"))
trv = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(trv)

TRAJ, LOG = "OUTPUT/trajectories-A.dat", "run_out.txt"
T_G, U_G, K_G, NU, D_P = trv.T_G, trv.U_G, trv.K_G, trv.NU, trv.D_P
TMIN, TMAX, T_0, H, slope = trv.TMIN, trv.TMAX, trv.T_0, trv.H, trv.slope
LZ, NZ = trv.LZ, trv.NZ

HALF_ULP, INT_FLOOR, T_BAND = 0.5e-6, 1.0e-8, 1.0   # F12.6 on x and T; integrator floor; resolved window
TOL_U, TOL_VW, TOL_RHO = 1.0e-4, 1.0e-6, 1.0e-6
MIN_PTS, N_PER_GROUP = 3, {210: 5, 300: 5, 750: 10, 840: 5}


def rho(T):
    """The table's density: the line 2950 - 0.5 (T - 250), which lookupTab continues outside the table."""
    return 2950.0 - 0.5 * (T - 250.0)


def d_of(m, T):
    return (6.0 * m / (math.pi * rho(T))) ** (1.0 / 3.0)


def half_ulp_e(v):
    """Half a unit in the last place of v printed E13.6 (0.dddddd x 10^e)."""
    return 0.5 * 10.0 ** (math.floor(math.log10(abs(v))) + 1 - 6)


def gauss_legendre(n):
    xs, ws = [], []
    for k in range(1, n + 1):
        x = math.cos(math.pi * (k - 0.25) / (n + 0.5))
        for _ in range(100):
            p0, p1 = 1.0, x
            for j in range(2, n + 1):
                p0, p1 = p1, ((2 * j - 1) * x * p1 - (j - 1) * p0) / j
            dp = n * (x * p1 - p0) / (x * x - 1.0)
            dx = p1 / dp
            x -= dx
            if abs(dx) < 1e-16:
                break
        xs.append(x)
        ws.append(2.0 / ((1.0 - x * x) * dp * dp))
    return xs, ws


GLX, GLW = gauss_legendre(20)


def dtdT(m, s, T):
    return s * m / (NU * K_G * math.pi * d_of(m, T) * abs(T_G - T))


def seg_slope(p, q):
    return slope(min(max(int(math.floor(0.5 * (p + q))), TMIN), TMAX - 1))


def piece(m, p, q):
    """Time to go from p to q inside one segment (either direction)."""
    s, a, b = seg_slope(p, q), min(p, q), max(p, q)
    h = 0.5 * (b - a)
    return h * sum(w * dtdT(m, s, 0.5 * (a + b) + h * x) for x, w in zip(GLX, GLW))


class Reference:
    """t(T) and its inverse from T0 toward T_g, stopping T_BAND short of T_g."""

    def __init__(self, T0, m):
        self.m, self.sign = m, (1.0 if T0 < T_G else -1.0)
        end = T_G - self.sign * T_BAND
        nodes = [T0]
        k = math.floor(T0) + 1 if self.sign > 0 else math.ceil(T0) - 1
        while (k - end) * self.sign < 0:
            if TMIN <= k <= TMAX and k != T0:
                nodes.append(float(k))
            k += int(self.sign)
        nodes.append(end)
        self.nodes, self.cum = nodes, [0.0]
        for p, q in zip(nodes[:-1], nodes[1:]):
            self.cum.append(self.cum[-1] + piece(m, p, q))

    def T_at(self, t):
        if t < 0.0 or t > self.cum[-1]:
            return None
        lo, hi = 0, len(self.cum) - 1
        while hi - lo > 1:
            mid = (lo + hi) // 2
            if self.cum[mid] <= t:
                lo = mid
            else:
                hi = mid
        a, b = self.nodes[lo], self.nodes[hi]
        base = self.cum[lo]
        for _ in range(200):                               # bisection inside the segment
            c = 0.5 * (a + b)
            if base + piece(self.m, self.nodes[lo], c) < t:
                a = c
            else:
                b = c
            if abs(b - a) < 1e-12:
                break
        return 0.5 * (a + b)

    def dTdx(self, T):
        s = slope(min(max(int(math.floor(T)), TMIN), TMAX - 1))
        return 1.0 / (U_G * dtdT(self.m, s, T))


def main():
    fails = []
    try:
        traj = trv.load_rows(TRAJ, 10)
        log = open(LOG, errors="replace").read()
        tmin, tmax, cols = load_properties(os.path.join(HERE, "INPUT", "properties.dat"))
    except (FileNotFoundError, ValueError, IndexError) as e:
        print(f"[FAIL] cannot read the outputs: {e}")
        return 1
    if (tmin, tmax) != (TMIN, TMAX) or max(abs(cols["rho"][T - TMIN] - rho(T)) for T in H) > TOL_RHO \
            or max(abs(cols["h"][T - TMIN] - H[T]) for T in H) > 1e-6:
        fails.append("G the fixture is not the recipe (range, density line or enthalpy column)")
    for w in ("Solving enthalpy equation", "Field mollification OFF"):
        if w not in log:
            fails.append(f"G witness: the log does not report '{w}'")

    refs, good = {}, {g: 0 for g in N_PER_GROUP}
    worst = {"V1": (0.0, ""), "D1": (0.0, ""), "D2": (0.0, ""), "M": (0.0, "")}

    def track(key, r, tol, what):
        if r / tol > worst[key][0]:
            worst[key] = (r / tol, what)
        if r > tol:
            fails.append(f"{key} {what}")

    for pid in sorted(traj):
        rows = sorted(traj[pid])
        x0, _, z0, _, _, _, Ta, d0p, m0p, _ = rows[0]
        kk = min(int(z0 / (LZ / NZ)) + 1, NZ)
        T0 = T_0[kk]
        group = int(round(T0))
        m = rho(T0) * math.pi / 6.0 * D_P ** 3
        if x0 != 0.0 or abs(Ta - T0) > HALF_ULP:
            fails.append(f"G ID={pid}: first row at x={x0}, T={Ta}, not x = 0, T = {T0}")
            continue
        ref = refs.setdefault(group, Reference(T0, m))
        n_win = 0
        for (x, y, z, u, v, w, T, d, mp, _) in rows:
            if abs(u - U_G) > TOL_U or abs(v) > TOL_VW or abs(w) > TOL_VW:
                fails.append(f"G ID={pid} x={x:.6f}: velocity ({u}, {v}, {w}) is not (u_g, 0, 0)")
                break
            track("M", abs(mp - m), half_ulp_e(m) + 1e-12 * m, f"ID={pid} x={x:.6f}: m={mp:.6e} vs {m:.6e}")
            ddT = d_of(m, T) * 0.5 / (3.0 * rho(T))
            track("D2", abs(d - d_of(m, T)), half_ulp_e(d) + ddT * HALF_ULP + 1e-12 * d,
                  f"ID={pid} x={x:.6f}: d={d:.6e} vs d(T)={d_of(m, T):.6e} at T={T:.6f}")
            if x <= 0.0 or abs(T - T_G) <= T_BAND:
                continue
            Tr = ref.T_at(x / U_G)
            if Tr is None:
                continue
            tolT = HALF_ULP * (1.0 + abs(ref.dTdx(T))) + INT_FLOOR
            n_win += 1
            track("V1", abs(T - Tr), tolT, f"ID={pid} x={x:.6f}: T={T:.6f} vs reference {Tr:.6f}")
            ddTr = d_of(m, Tr) * 0.5 / (3.0 * rho(Tr))
            track("D1", abs(d - d_of(m, Tr)), half_ulp_e(d) + ddTr * tolT + 1e-12 * d,
                  f"ID={pid} x={x:.6f}: d={d:.6e} vs d(T_ref)={d_of(m, Tr):.6e}")
        if n_win >= MIN_PTS:
            good[group] += 1

    for g, n in N_PER_GROUP.items():
        if good[g] < n:
            fails.append(f"G non-vacuity: {good[g]} parcels injected at {g} K with >= {MIN_PTS} window rows (need {n})")
    for k in ("V1", "D1", "D2", "M"):
        print(f"{k} worst: {worst[k][1]} (resid/tol {worst[k][0]:.3f})")
    print("verified parcels per group: " + ", ".join(f"{g} K {good[g]}/{n}" for g, n in N_PER_GROUP.items()))
    for f in fails[:40]:
        print(f"[FAIL] {f}")
    if fails:
        print(f"\n[FAIL] {len(fails)} violation(s)")
        return 1
    print("\n[PASS] with a varying density the printed diameter follows d(T) at constant mass, and the temperature "
          "follows the energy balance with that diameter.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
