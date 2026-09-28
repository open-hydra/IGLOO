#!/usr/bin/env python3
"""The property table INPUT/<prefix>properties.dat as IGLOO reads it, for oracles.

    from proptab import load_properties, lookup, table_value
    Tmin, Tmax, cols = load_properties("INPUT/properties.dat")   # zone 1
    h = lookup(cols["h"], Tmin, Tmax, Tp)                        # lookupTab: h, sigma, mu
    rho = table_value(cols["rho"], Tmin, Tmax, Tp)               # tableValue: rho, psat

The header is the first line holding VARIABLES; its double-quoted names name the columns, the
first being the temperature. Known names map to the slots IGLOO reads (T, cp, rho, h, psat) with
the alias set of the solver's reader; any other column keeps its own name. Rows are consecutive
integer kelvins, one zone per material in phase-file order.

The two interpolations replicate the solver statement for statement: `lookup` clamps the INDEX to
[Tmin, Tmax-1] and so extrapolates linearly outside the table; `table_value` clamps the
TEMPERATURE and so takes the end value outside it. An oracle must use the one the solver uses for
the property it checks.
"""
import math

SLOTS = {
    "Temperature": "T", "temperature": "T",
    "Cp": "cp", "cp": "cp",
    "Density": "rho", "density": "rho",
    "Enthalpy": "h", "enthalpy": "h", "Enthalpy_abs": "h", "enthalpy_abs": "h",
    "Psat": "psat", "psat": "psat", "PSAT": "psat",
}


def header_names(path):
    """The quoted names of the first VARIABLES line, in order."""
    with open(path) as f:
        for line in f:
            if "VARIABLES" in line:
                return line.split('"')[1::2]
    raise ValueError(f"{path}: no VARIABLES line")


def load_properties(path, zone=1):
    """(Tmin, Tmax, {slot or name: [values]}) of one zone; also 'datum' = relative|absolute."""
    names = header_names(path)
    keys = [SLOTS.get(n, n) for n in names]
    zones, cur = [], None
    with open(path) as f:
        for line in f:
            tok = line.split()
            try:
                vals = [float(t) for t in tok[:len(names)]]
            except ValueError:
                vals = None
            if tok and vals is not None and len(vals) == len(names):
                if cur is None:
                    cur = []
                    zones.append(cur)
                cur.append(vals)
            else:
                cur = None
    rows = zones[zone - 1]
    cols = {k: [r[i] for r in rows] for i, k in enumerate(keys)}
    Tmin = int(round(cols["T"][0]))
    Tmax = Tmin + len(rows) - 1
    cols["datum"] = "absolute" if any(n.lower() == "enthalpy_abs" for n in names) else "relative"
    return Tmin, Tmax, cols


def lookup(tab, Tmin, Tmax, T):
    """lookupTab (enthalpy, sigma, mu): the node index clamped to [Tmin, Tmax-1], linear in T from it."""
    Ti = min(max(int(T), Tmin), Tmax - 1)
    v0 = tab[Ti - Tmin]
    v1 = tab[Ti + 1 - Tmin]
    return v0 + (v1 - v0) * (T - Ti)


def table_value(tab, Tmin, Tmax, T):
    """tableValue (density, psat): linear between the nodes, the end value outside, NaN for a NaN temperature."""
    if math.isnan(T):
        return T
    Tc = min(max(T, float(Tmin)), float(Tmax))
    if Tc >= Tmax:
        return tab[Tmax - Tmin]
    i = int(Tc)
    return tab[i - Tmin] + (tab[i + 1 - Tmin] - tab[i - Tmin]) * (Tc - i)
