#!/usr/bin/env python3
"""Read an IGLOO .dat file (trajectories, exit records, scatter, snapshots) by the names of its header.

zones(path) returns one dict per Zone -- {"title": <Zone T= text>, "names": [<variables= names>],
"rows": [{name: value}]} -- an integer token (the ID, a boundary code) read as an int, every other as a
float. A row whose width
differs from its header's name count raises ValueError: a checker that reads columns by name sees a layout
change in its header, never by guessing a width. rows(path) concatenates the rows of every zone.
solid_phase(f, T, Tmelt) is the phase of a solidifying material's row: the files carry f, not the phase.

>>> import os, tempfile
>>> d = tempfile.mkdtemp()
>>> p = os.path.join(d, "outloc-A.dat")
>>> _ = open(p, "w").write(
...     ' variables="X","Y","Z","T","|u<sub>p</sub>|","<greek>a</greek>","mdot","Af","ID"\\n'
...     'Zone T="Mat A Group  1"\\n'
...     '    0.150000    0.005000    0.035000  600.000000    9.999968    0.000001 0.618182E-03 0.100000E-03      16\\n'
...     'Zone T="Mat A Group  2"\\n'
...     '    0.150000    0.045000    0.005000  600.000000    9.999968    0.000001 0.618182E-03 0.100000E-03      21\\n')
>>> z = zones(p)
>>> [x["title"] for x in z]
['Mat A Group  1', 'Mat A Group  2']
>>> z[0]["names"][-2:]
['Af', 'ID']
>>> r = z[0]["rows"][0]
>>> r["ID"], r["mdot"], r["X"]
(16, 0.000618182, 0.15)
>>> [r["ID"] for r in rows(p)]
[16, 21]
>>> _ = open(p, "a").write('    0.150000    0.045000      22\\n')
>>> zones(p)
Traceback (most recent call last):
    ...
ValueError: outloc-A.dat:6: 3 values for 9 names
>>> _ = open(p, "w").write(
...     ' variables="X","Y","Z","T","|u<sub>p</sub>|","<greek>a</greek>","mdot","Af","ID","t","bc","d<sub>p</sub>"\\n'
...     'Zone T="Mat A Group  1"\\n'
...     '    0.150000    0.005000    0.035000  600.000000    9.999968    0.000001 0.618182E-03 0.100000E-03      16'
...     '  1.50000000E-02  400 0.118900E-04\\n')
>>> r = rows(p)[0]
>>> r["t"], r["bc"], r["d<sub>p</sub>"], r["ID"]
(0.015, 400, 1.189e-05, 16)
"""
import os
import re

INT = re.compile(r"[+-]?\d+$")


def _names(line):
    return re.findall(r'"([^"]*)"', line)


def zones(path):
    """One dict per Zone: title, the names of the variables= header in force, and the rows keyed by them."""
    out, names, cur = [], None, None
    with open(path) as fh:
        for n, line in enumerate(fh, 1):
            s = line.strip()
            if not s:
                continue
            low = s.lower()
            if low.startswith("variables"):
                names = _names(s)
                continue
            if low.startswith("zone"):
                m = re.search(r'T\s*=\s*"([^"]*)"', s)
                cur = {"title": m.group(1) if m else "", "names": list(names or []), "rows": []}
                out.append(cur)
                continue
            if low.startswith("title"):
                continue
            if names is None:
                raise ValueError("%s:%d: a data row before the variables= header" % (os.path.basename(path), n))
            if cur is None:
                cur = {"title": "", "names": list(names), "rows": []}
                out.append(cur)
            vals = s.split()
            if len(vals) != len(cur["names"]):
                raise ValueError("%s:%d: %d values for %d names" % (os.path.basename(path), n, len(vals),
                                                                    len(cur["names"])))
            cur["rows"].append({k: (int(v) if INT.match(v) else float(v)) for k, v in zip(cur["names"], vals)})
    return out


def rows(path):
    """The rows of every zone of the file, in file order."""
    return [r for z in zones(path) for r in z["rows"]]


def solid_phase(f, T, Tmelt):
    """The phase of a solidifying material's row from its frozen fraction f and temperature T: 3 solid
    (f = 1), 2 plateau (0 < f < 1), else 1 undercooled liquid (T < Tmelt) or 0 liquid (T >= Tmelt).

    >>> [solid_phase(f, T, 2327.0) for f, T in ((0.0, 2400.0), (0.0, 2327.0), (0.0, 2000.0),
    ...                                         (0.54369, 2327.0), (1.0, 2327.0), (1.0, 907.3))]
    [0, 0, 1, 2, 3, 3]
    """
    if f >= 1.0:
        return 3
    if f > 0.0:
        return 2
    return 1 if T < Tmelt else 0


if __name__ == "__main__":
    import doctest
    import sys
    failed, tried = doctest.testmod()
    print("datcols: %d doctest examples, %d failed" % (tried, failed))
    sys.exit(1 if failed else 0)
