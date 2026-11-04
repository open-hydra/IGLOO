#!/usr/bin/env python3
"""The enthalpy state of a parcel along its measured path, for the oracles of a case whose particle state is the
enthalpy (a tabulated, varying cp).

    tab = Table("INPUT/properties.dat")   # h(T) as lookupTab, T(h) as comp_TfromTab, cp(T) the slope of h(T)
    ys = march(ts, drivers, y0, rate, step, nsub, reset)

`march` integrates the state y (a list) by RK4 from y0 at ts[0] and returns {row index: y} at every step-th row and
the last. rate(t, y, v) is dy/dt, v each driver ({name: [value per row]}) at t: the cubic Lagrange polynomial through
the four nearest rows used (fewer when the rows are fewer). reset(k, y) returns the state an interval starts from at
row k (a component re-anchored on the row's measured value); nsub(dt, t, y, v) is the RK4 step count of the interval
(at least 8 per row it spans). `runs(values)` splits a row list where a driver changes between consecutive rows.
"""
import math

from proptab import load_properties, lookup


class Table:
    """The h column of a property table, as the solver interpolates and inverts it."""

    def __init__(self, path):
        self.tmin, self.tmax, cols = load_properties(path)
        self.hcol = cols["h"]
        self.cols = cols

    def h(self, T):
        return lookup(self.hcol, self.tmin, self.tmax, T)

    def T(self, h):
        """comp_TfromTab: the end segments extended, the bracketing node by bisection."""
        tab, lo, hi = self.hcol, self.tmin, self.tmax
        if h < tab[0]:
            return lo + (h - tab[0]) / (tab[1] - tab[0])
        if h >= tab[-1]:
            return hi + (h - tab[-1]) / (tab[-1] - tab[-2])
        i, j = lo, hi
        while j - i > 1:
            k = (i + j) // 2
            if tab[k - lo] <= h:
                i = k
            else:
                j = k
        return i + (h - tab[i - lo]) / (tab[i + 1 - lo] - tab[i - lo])

    def cp(self, T):
        """The slope of h(T) on the table segment that holds T (lookupTab's index clamp)."""
        i = min(max(int(math.floor(T)), self.tmin), self.tmax - 1)
        return self.hcol[i + 1 - self.tmin] - self.hcol[i - self.tmin]


def runs(values):
    """[(first, last)] row ranges over which `values` does not change."""
    out, start = [], 0
    for k in range(1, len(values)):
        if values[k] != values[k - 1]:
            out.append((start, k - 1))
            start = k
    out.append((start, len(values) - 1))
    return out


def _lagrange(tu, yu, nodes, t):
    q = 0.0
    for p in nodes:
        w = 1.0
        for r in nodes:
            if r != p:
                w *= (t - tu[r]) / (tu[p] - tu[r])
        q += w * yu[p]
    return q


def march(ts, drivers, y0, rate, step=1, nsub=None, reset=None):
    """{row index: y} at every step-th row and the last, from y0 at row 0 (see the module docstring)."""
    idx = list(range(0, len(ts), step))
    if idx[-1] != len(ts) - 1:
        idx.append(len(ts) - 1)
    tu = [ts[i] for i in idx]
    du = {k: [v[i] for i in idx] for k, v in drivers.items()}
    npt = len(idx)
    y = list(y0)
    out = {0: list(y)}
    for j, (a, b) in enumerate(zip(idx[:-1], idx[1:])):
        if npt >= 4:
            lo = min(max(j - 1, 0), npt - 4)
            nodes = range(lo, lo + 4)
        else:
            nodes = range(npt)

        def vals(t, nodes=nodes):
            return {k: _lagrange(tu, du[k], nodes, t) for k in du}

        if reset is not None:
            y = reset(a, y)
        ta, tb = ts[a], ts[b]
        n = 8 * (b - a)
        if nsub is not None:
            n = max(n, nsub(tb - ta, ta, y, vals(ta)))
        dt = (tb - ta) / n
        t = ta
        for _ in range(n):
            v1, v2, v4 = vals(t), vals(t + 0.5 * dt), vals(t + dt)
            k1 = rate(t, y, v1)
            k2 = rate(t + 0.5 * dt, [p + 0.5 * dt * q for p, q in zip(y, k1)], v2)
            k3 = rate(t + 0.5 * dt, [p + 0.5 * dt * q for p, q in zip(y, k2)], v2)
            k4 = rate(t + dt, [p + dt * q for p, q in zip(y, k3)], v4)
            y = [p + dt / 6.0 * (q1 + 2.0 * q2 + 2.0 * q3 + q4) for p, q1, q2, q3, q4 in zip(y, k1, k2, k3, k4)]
            t += dt
        out[b] = list(y)
    return out
