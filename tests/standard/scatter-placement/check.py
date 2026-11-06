#!/usr/bin/env python3
"""scatter-placement: where the scatter cloud's markers fall, how many there are and what they weigh, on khrt-stress's
box (KH-RT breakup with mShedLim 0.01: hundreds of children, Rayleigh-Taylor shatters that multiply a parcel's number
rate), out-file ALL.

A parcel writes markers at a constant rate r = n_ref/dNscat, n_ref the number rate it entered with raised to the
material's floor (a breakup child keeps its parent's), at the instants its accumulator crosses a whole number, inside
the solver's steps: the k-th marker at t_k = t_birth + (k - frac(ID phi))/r. Each marker carries its weight
w = np(t_k)/r (column "weight"), the droplets it stands for, np interpolated linearly in its step. The floor is
F nbar, F the [IGLOO-General] scatter-floor (0.1 when absent) and nbar = Ndot/nStreams the mean of the 25 injected
streams' rates; a rate below (1 - 1e-6) F nbar takes the floor.

  P1  markers within TOL = 1e-6 m of an x-face plane (x = i dx) are at most 3 x 2 TOL/dx of all markers: a marker
      uniform in time is near a face as often as the face's share of the cell (0.0008 here), while a marker at a
      step end crowds the faces, where the steps are cut
  P2  no two markers alike (every printed column): a marker at the solver's restart point repeats the step end
  P3  the markers' fractional x in their cells against U(0,1): Kolmogorov-Smirnov p >= 1e-3
  C   per parcel, the marker count N = floor(frac(ID phi) + r tau), tau its exit t less its birth t, r with the floor
  F   the floor's premise: nbar/dNscat = fsample-traj V_ref/L_ref (the quantum's definition, nbar taken over the 25
      roots), V_ref the roots' mean first-row speed and L_ref the mesh's bounding-box diagonal
  W1  the Euler moment, the sum over cells of n_p V (euler1.tec, the ODE's integral of np dt), lies in the bracket
      the parcels' printed states give [sum I_lo, sum I_hi]
  W2  per parcel, the markers' sum of w lies in [I_lo - B, I_hi + B], B = (TV + M)/r
  W3  over all parcels, the sum of w differs from the Euler moment by at most the sum of the parcels' bounds on
      |S_p - I_p|: I_hi for a parcel without a marker, B + I_hi - I_lo for the others

C: the accumulator starts at frac(ID phi) and adds r h on every step the parcel keeps, and every crossing of a whole
number writes a marker, the exit face's included, so N = floor(frac(ID phi) + r tau) whatever the steps. r is the
number rate of the parcel's injected ancestor at its entry, np = mdot/m from its birth flow and its first row, raised
to the floor F nbar when below (1 - 1e-6) of it, over dNscat (the children are all of the first generation, and the
case has no residue); nbar is the mean of the 25 roots' rates, an interval from theirs, and a rate whose interval
meets the floor's accepts either branch. Every printed input enters with
its half-ULP and a parcel whose interval of frac(ID phi) + r tau holds an integer accepts either floor (dNscat's
ES10.3 dominates the interval); at least C_ONE of the parcels have a window of one value, where N - 1 and N + 1 both
fail. The weight clauses gate the np factor of w, C gates r: a wrong rate leaves the sum of w a Riemann sum of the
same integral, a wrong np leaves the count.

W: the bracket and the bound come from the model and the parcel's printed states (its trajectory rows, out-time on,
and its exit record), not from the markers. With a constant density and no evaporation np = mdot/m, the flow mdot
falls only at a shed (by the child's birth flow: a parcel's birth flow is its exit flow plus its children's) and the
droplet mass m never grows (KH stripping shrinks it at a constant flow, a shatter divides it at a constant flow, a
shed keeps it); W checks that premise on the printed m. On the interval between two printed states np therefore lies
between mdot(end)/m(start) and mdot(start)/m(end): I_lo and I_hi are those values times the interval, and the
interval's variation of np is at most mdot(start)(1/m(end) - 1/m(start)) plus its sheds over m(end). TV is their
sum over the parcel's life, M the largest upper value. The markers are the Riemann sum, step 1/r, of the linear
interpolant of np over the solver's steps, whose integral lies in [I_lo, I_hi] and whose variation and maximum are
at most TV and M; the sum's error is at most (TV + M)/r (TV for the samples inside their intervals, a shatter's
jump included; M for the seed and the unwritten tail). The printed inputs' rounding enters as factors (dNscat's
ES10.3 dominates B) and the weights with their half-ULP. W3 is loose where parcels hold few markers (their B exceeds
their integral); W2 holds each parcel to its own.

Non-vacuous: P1-P3 fail on markers at step ends (the solver's restart points and the steps cut at the faces); C fails
on a child that takes its own number rate at the shed, or its parent's np at the shed, instead of its lineage's
rate, and on a root or a child whose rate misses the floor; W2 and W3 fail on weights frozen at the parcel's entry
(np's growth through the shatters escapes them) and on weights left at the entry rate under a floored marker rate.
"""
import math
import os
import re
import sys
from collections import Counter, defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "tools"))
import datcols                                  # noqa: E402
import ks_ad                                    # noqa: E402

RHO = 1000.0                 # [GPB-Phase1] rho
DX = 0.0025                  # the box's cell along x (solfile.tec)
TOL = 1.0e-6                 # the F12.6 print resolution [m]
FACE_MAX = 3.0 * 2.0 * TOL / DX
P_MIN = 1.0e-3
PHI = 0.6180339887498949     # the accumulator's per-ID phase (Lib_Integration.f90, integrate)
C_ONE = 0.9                  # C's share of one-valued windows (measured 721 of 726)
R_I = 1.0e-4                 # the bracket's printed inputs: m, mdot (E13.6E2), t (ES16.8E2)
R_B = 2.0e-3                 # B's: dNscat (ES10.3, 5e-4), m through TV (5e-6 per shed piece)
R_M = 1.0e-5                 # m between rows (E13.6E2)
R_MX = 2.0e-5                # m at the exit record, from d (E13.6E2, cubed)
FLOOR_DEFAULT = 0.1          # [IGLOO-General] scatter-floor when absent
FLOOR_MARGIN = 1.0e-6        # a rate takes the floor only below (1 - FLOOR_MARGIN) of it
TRAJ_SAMPLE = 100            # [IGLOO-General] fsample-traj (absent: 100)


def fail(msg):
    print("[FAIL] " + msg)
    return 1


def half_ulp_e(v):
    """Half-ULP of an E13.6 value (0.dddddd E+ee)."""
    return 0.0 if v == 0.0 else 0.5e-6 * 10.0 ** (math.floor(math.log10(abs(v))) + 1)


def half_ulp_es(text):
    """Half-ULP of a printed ES value, from its own digits; a printed zero is exact."""
    m = re.match(r"\s*[-+]?\d\.(\d+)[Ee]([-+]?\d+)", text)
    return 0.0 if float(text) == 0.0 else 0.5 * 10.0 ** (int(m.group(2)) - len(m.group(1)))


def es8(v):
    """The ES16.8E2 text of a printed t (a 9-digit decimal: the float holds it exactly)."""
    return "%.8E" % v


def cell_volumes(path):
    """Sum of n_p V over the cells of every block of a BLOCK-packed euler<fam>.tec (nodal x y z, cell-centred n_p)."""
    txt = open(path).read().split("\n")
    names = txt[0].split("=", 1)[1].split()
    ip = names.index("n<sub>p")
    zones, cur = [], None
    for ln in txt[1:]:
        s = ln.strip()
        if s.upper().startswith("ZONE"):
            m = re.search(r"I=(\d+), J=(\d+), K=(\d+)", s)
            cur = [int(m.group(1)), int(m.group(2)), int(m.group(3)), []]
            zones.append(cur)
        elif s:
            cur[3].append(float(s))
    total = 0.0
    for ni, nj, nk, v in zones:
        nn, nc = ni * nj * nk, (ni - 1) * (nj - 1) * (nk - 1)
        x, y, z = v[0:nn], v[nn:2 * nn], v[2 * nn:3 * nn]
        npv = v[3 * nn + (ip - 3) * nc: 3 * nn + (ip - 2) * nc]

        def node(i, j, k):
            return i + ni * (j + nj * k)
        c = 0
        for k in range(nk - 1):
            for j in range(nj - 1):
                for i in range(ni - 1):
                    a, b = node(i, j, k), node(i + 1, j + 1, k + 1)
                    vol = (x[b] - x[a]) * (y[b] - y[a]) * (z[b] - z[a])    # the box's cells are axis-aligned
                    total += npv[c] * vol
                    c += 1
    return total


def scatter_floor():
    """[IGLOO-General] scatter-floor of input.ini, the default when absent."""
    sec = None
    for ln in open("input.ini"):
        s = ln.split(";")[0].strip()
        if s.startswith("["):
            sec = s
        elif sec == "[IGLOO-General]" and "=" in s and s.split("=", 1)[0].strip() == "scatter-floor":
            return float(s.split("=", 1)[1])
    return FLOOR_DEFAULT


def mesh_diagonal(path):
    """Bounding-box diagonal of the nodes of every zone of a BLOCK-packed gas file (x, y, z first)."""
    lo, hi = [math.inf] * 3, [-math.inf] * 3
    txt = open(path).read().split("\n")
    zones = [i for i, ln in enumerate(txt) if ln.strip().startswith("Zone")]
    for zi, i0 in enumerate(zones):
        m = re.search(r"I=(\d+),\s*J=(\d+),\s*K=(\d+)", txt[i0])
        nn = int(m.group(1)) * int(m.group(2)) * int(m.group(3))
        i1 = zones[zi + 1] if zi + 1 < len(zones) else len(txt)
        vals = []
        for ln in txt[i0 + 1:i1]:
            if "=" in ln:
                continue
            vals.extend(float(v) for v in ln.split())
            if len(vals) >= 3 * nn:
                break
        for c in range(3):
            lo[c] = min(lo[c], min(vals[c * nn:(c + 1) * nn]))
            hi[c] = max(hi[c], max(vals[c * nn:(c + 1) * nn]))
    return math.sqrt(sum((h - l) ** 2 for l, h in zip(lo, hi)))


def lineage():
    """The parcels' printed states and breakup tree: dNscat, the exit records, the trajectory knots (t, m), each
    parcel's parent, children and injected root, and its birth flow (its exit flow plus its children's) with the
    sum of the printed flows' half-ULPs; a string on a malformed tree."""
    m = re.search(r"scatter cloud \[A\]: dNscat=\s*(\S+)", open("run_out.txt").read())
    exits = {r["ID"]: r for r in datcols.rows("OUTPUT/outloc-A.dat")}
    rows, parent, v0 = {}, {}, {}
    for z in datcols.parcel_zones("OUTPUT/trajectories-A.dat"):
        rows.setdefault(z["id"], []).extend((r["t"], r["m<sub>p"]) for r in z["rows"])
        parent[z["id"]] = z["parent"]
        v0.setdefault(z["id"], tuple(z["rows"][0][c] for c in ("U", "V", "W")))
    ids = sorted(exits)
    if set(ids) != set(rows):
        return "%d exit records, %d parcels with trajectory rows" % (len(ids), len(rows))
    if any(parent[i] >= i for i in ids if parent[i]):
        return "a child's ID is not above its parent's"
    kids = defaultdict(list)
    for i in ids:
        if parent[i]:
            kids[parent[i]].append(i)
    mdot_b, mdot_h, root = {}, {}, {}
    for i in sorted(ids, reverse=True):
        mdot_b[i] = exits[i]["mdot"] + math.fsum(mdot_b[c] for c in kids[i])
        mdot_h[i] = half_ulp_e(exits[i]["mdot"]) + math.fsum(mdot_h[c] for c in kids[i])
    for i in ids:
        root[i] = root[parent[i]] if parent[i] else i
        if root[i] == i and rows[i][0][0] != 0.0:
            return "ID %d has no parent and is born at t = %s, not injected" % (i, es8(rows[i][0][0]))
    return {"dn": float(m.group(1)), "dn_h": half_ulp_es(m.group(1)), "exits": exits, "rows": rows, "ids": ids,
            "parent": parent, "kids": kids, "root": root, "mdot_b": mdot_b, "mdot_h": mdot_h, "v0": v0}


def rates(lin, floor):
    """Each root's marker number rate with the floor: {root: (lo, hi, central)}, and the floor's (lo, hi, central).
    A root's own rate np = mdot/m from its birth flow and first row (each with its half-ULP); nbar the mean over the
    roots; a rate wholly below (1 - FLOOR_MARGIN) of the floor's interval takes the floor, one wholly at or above it
    keeps its own, one that meets it accepts either."""
    own = {}
    for q in lin["ids"]:
        if lin["root"][q] != q:
            continue
        m0, f0, fh = lin["rows"][q][0][1], lin["mdot_b"][q], lin["mdot_h"][q]
        own[q] = ((f0 - fh) / (m0 + half_ulp_e(m0)), (f0 + fh) / (m0 - half_ulp_e(m0)), f0 / m0)
    n = len(own)
    fl = tuple(floor * math.fsum(r[k] for r in own.values()) / n for k in range(3))
    out = {}
    for q, (lo, hi, c) in own.items():
        if hi < (1.0 - FLOOR_MARGIN) * fl[0]:
            out[q] = fl
        elif lo >= (1.0 - FLOOR_MARGIN) * fl[1]:
            out[q] = (lo, hi, c)
        else:
            out[q] = (min(lo, fl[0]), max(hi, fl[1]), fl[2] if c < (1.0 - FLOOR_MARGIN) * fl[2] else c)
    return out, fl, own


def premise(lin, own, floor):
    """F: nbar/dNscat = fsample-traj V_ref/L_ref, nbar the roots' mean rate, V_ref their mean first-row speed."""
    n = len(own)
    nbar = math.fsum(r[2] for r in own.values()) / n
    nbar_h = math.fsum((r[1] - r[0]) / 2.0 for r in own.values()) / n
    v = [math.sqrt(sum(c * c for c in lin["v0"][q])) for q in own]
    v_ref = math.fsum(v) / n
    v_h = 0.5e-6 * math.sqrt(3.0)                   # F12.6 velocity components
    l_ref = mesh_diagonal(os.path.join("INPUT", "solfile.tec"))
    got, want = nbar / lin["dn"], TRAJ_SAMPLE * v_ref / l_ref
    tol = got * (lin["dn_h"] / lin["dn"] + nbar_h / nbar) + want * (v_h / v_ref + 1.0e-6)
    msg = "F nbar/dNscat %.6e against fsample-traj V_ref/L_ref %.6e (nbar %.6e 1/s over %d roots, V_ref %.6f m/s, " \
          "L_ref %.6f m, floor %g nbar), tol %.2e" % (got, want, nbar, n, v_ref, l_ref, floor, tol)
    if abs(got - want) <= tol:
        print(msg + " PASS")
        return 0
    return fail(msg)


def placement(scat):
    rc, n = 0, len(scat)
    near, frac = 0, []
    for r in scat:
        u = r["X"] / DX
        f = u - math.floor(u)
        frac.append(f)
        if min(f, 1.0 - f) * DX <= TOL:
            near += 1
    share = near / n
    msg = "P1 %d of %d markers (%.4f) within %.0e m of an x-face, bound %.4f" % (near, n, share, TOL, FACE_MAX)
    if share > FACE_MAX:
        rc |= fail(msg)
    else:
        print(msg + " PASS")
    dup = n - len(set(tuple(sorted(r.items())) for r in scat))
    if dup:
        rc |= fail("P2 %d of %d markers repeat another (%.4f)" % (dup, n, dup / n))
    else:
        print("P2 no two of the %d markers alike PASS" % n)
    d = ks_ad.ks_d(frac, lambda x: min(max(x, 0.0), 1.0))
    p = ks_ad.ks_pvalue(d, n)
    msg = "P3 fractional x against U(0,1): D = %.4f over %d markers, p = %.3g" % (d, n, p)
    if p < P_MIN:
        rc |= fail(msg + " < %.0e" % P_MIN)
    else:
        print(msg + " PASS")
    return rc


def counts(scat, lin, rate):
    dn, dh, exits, rows, ids = lin["dn"], lin["dn_h"], lin["exits"], lin["rows"], lin["ids"]
    n = Counter(r["ID"] for r in scat)
    out, two = [], []
    d_lo, d_hi = dn - dh, dn + dh
    for i in ids:
        np_lo, np_hi, _ = rate[lin["root"][i]]
        te, tb = exits[i]["t"], rows[i][0][0]
        th = half_ulp_es(es8(te)) + half_ulp_es(es8(tb))
        r_lo, r_hi = np_lo * max(te - tb - th, 0.0), np_hi * (te - tb + th)
        ph = math.fmod(i * PHI, 1.0)
        a, b = math.floor(ph + r_lo / (dn + dh)), math.floor(ph + r_hi / (dn - dh))
        if a != b:
            two.append(i)
        if not a <= n[i] <= b:
            out.append((i, n[i], a, b))
        d_lo = max(d_lo, r_lo / (n[i] + 1 - ph))
        if n[i] > ph:
            d_hi = min(d_hi, r_hi / (n[i] - ph))
    for i, k, a, b in out[:10]:
        print("   ID %d (parent %d): %d markers, floor(frac(ID phi) + r tau) in [%d, %d]" % (
            i, lin["parent"][i], k, a, b))
    msg = "C %d of %d parcels off floor(frac(ID phi) + r tau) (%d markers); %d windows of one value, two on IDs %s" % (
        len(out), len(ids), len(scat), len(ids) - len(two), " ".join(str(i) for i in two[:12]) or "none")
    if out:
        rc = fail(msg)
    elif len(ids) - len(two) < C_ONE * len(ids):
        rc = fail(msg + ": fewer than %.0f %% of one value" % (100.0 * C_ONE))
    else:
        rc = 0
        print(msg + " PASS")
    if d_lo < d_hi:
        print("   one quantum fits every count: dNscat in (%.7g, %.7g], printed %.4g +- %.0e" % (d_lo, d_hi, dn, dh))
    else:
        print("   no quantum in the printed dNscat %.4g +- %.0e fits every count" % (dn, dh))
    return rc


def weights(scat, lin, rate):
    rc = 0
    if "weight" not in scat[0]:
        return fail("W the scatter file has no weight column")
    dn, exits, rows, ids, kids, mdot_b = lin["dn"], lin["exits"], lin["rows"], lin["ids"], lin["kids"], lin["mdot_b"]
    s_p, h_p = defaultdict(float), defaultdict(float)
    for r in scat:
        s_p[r["ID"]] += r["weight"]
        h_p[r["ID"]] += half_ulp_e(r["weight"])
    lo_t = hi_t = g_b = 0.0
    grows, out, tight = 0, [], [0, 0.0]
    for i in ids:
        knots = rows[i] + [(exits[i]["t"], RHO * math.pi * exits[i]["d<sub>p</sub>"] ** 3 / 6.0)]
        sheds = sorted((rows[c][0][0], mdot_b[c]) for c in kids[i])

        def after(t):                                # mdot just after t: the exit flow and the sheds still to come
            return exits[i]["mdot"] + math.fsum(f for ts, f in sheds if ts > t)
        lo = hi = tv = mx = 0.0
        for j, ((ta, ma), (tb, mb)) in enumerate(zip(knots, knots[1:])):
            if mb > ma * (1.0 + (R_MX if j == len(knots) - 2 else R_M)):
                grows += 1
            fa, fb, h = after(ta), after(tb), tb - ta
            lo += fb / ma * h
            hi += fa / mb * h
            tv += fa * max(0.0, 1.0 / mb - 1.0 / ma) + math.fsum(f for ts, f in sheds if ta < ts <= tb) / mb
            mx = max(mx, fa / mb)
        b = dn * (tv + mx) / rate[lin["root"][i]][2] * (1.0 + R_B)
        a, c = lo * (1.0 - R_I) - b - h_p[i], hi * (1.0 + R_I) + b + h_p[i]
        if not a <= s_p[i] <= c:
            out.append((i, s_p[i], lo, hi, b))
        if c - a <= lo:
            tight[0] += 1
            tight[1] += lo
        lo_t += lo
        hi_t += hi
        g_b += hi * (1.0 + R_I) if s_p[i] == 0.0 else b + hi * (1.0 + R_I) - lo * (1.0 - R_I) + h_p[i]
    if grows:
        rc |= fail("W the bound's premise: the droplet mass grows on %d intervals between printed states" % grows)
    e = cell_volumes("OUTPUT/euler1.tec")
    msg = "W1 sum of n_p V %.7e in the parcels' bracket [%.7e, %.7e] (width %.3e of the moment)" % (
        e, lo_t, hi_t, (hi_t - lo_t) / e)
    if lo_t * (1.0 - R_I) <= e <= hi_t * (1.0 + R_I):
        print(msg + " PASS")
    else:
        rc |= fail(msg)
    for i, s, lo, hi, b in out[:10]:
        print("   ID %d: sum of w %.6e, bracket [%.6e, %.6e], bound %.4e" % (i, s, lo, hi, b))
    msg = "W2 %d of %d parcels outside [I_lo - B, I_hi + B]; on %d (%.3f of the moment) the window is within +-50 %%" \
        % (len(out), len(ids), tight[0], tight[1] / lo_t)
    if out:
        rc |= fail(msg)
    else:
        print(msg + " PASS")
    s = math.fsum(s_p.values())
    msg = "W3 sum of w %.6e against the moment: %.3e of it, bound %.3e of it" % (s, (s - e) / e, g_b / e)
    if abs(s - e) <= g_b:
        print(msg + " PASS")
    else:
        rc |= fail(msg)
    return rc


def main():
    scat = datcols.rows("OUTPUT/scatter-A.dat")
    if not scat:
        return fail("no scatter rows")
    rc = placement(scat)
    lin = lineage()
    if isinstance(lin, str):
        rc |= fail("the breakup tree: " + lin)
    else:
        floor = scatter_floor()
        rate, fl, own = rates(lin, floor)
        nfl = sum(1 for q in own if own[q][1] < (1.0 - FLOOR_MARGIN) * fl[0])
        print("   floor %g nbar = %.6e 1/s: %d of %d roots below it" % (floor, fl[2], nfl, len(own)))
        rc |= premise(lin, own, floor)
        rc |= counts(scat, lin, rate)
        rc |= weights(scat, lin, rate)
    print("PASS" if rc == 0 else "FAIL")
    return rc


if __name__ == "__main__":
    sys.exit(main())
