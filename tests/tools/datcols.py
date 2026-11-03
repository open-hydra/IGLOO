#!/usr/bin/env python3
"""Read an IGLOO .dat file (trajectories, exit records, scatter, snapshots) by the names of its header.

zones(path) returns one dict per Zone -- {"title": <Zone T= text>, "names": [<variables= names>],
"rows": [{name: value}], "I": <the zone line's I=, None without one>} -- an integer token (the ID, a boundary
code) read as an int, every other as a float. A row whose width
differs from its header's name count raises ValueError: a checker that reads columns by name sees a layout
change in its header, never by guessing a width. rows(path) concatenates the rows of every zone.
group(title) is the group of a zone title. A trajectory file holds one zone per parcel image, titled
`Mat <mat> Group <g> ID <id> parent <pid>[ image <k>]`: parcel(title) reads those fields, and
parcel_zones(path) returns a trajectory file's zones with them, after checking the file's form.
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
>>> group(z[1]["title"]), group("Mat A Group100 ID 3 parent 0")
(2, 100)
>>> parcel("Mat A Group  1 ID 7 parent 2 image 3"), parcel("Mat A Group  1 ID 7 parent 0"), parcel("Mat A Group  1")
(('A', 1, 7, 2, 3), ('A', 1, 7, 0, 0), None)
>>> t = os.path.join(d, "trajectories-A.dat")
>>> head = ' variables="X","Y","Z","U","V","W","T","d<sub>p","m<sub>p","ID"\\n'
>>> row = '    0.000000    0.005000    0.005000   10.000000    0.000000    0.000000  600.000000 0.1E-04 0.2E-11 %7d\\n'
>>> zl = 'ZONE T="Mat A Group  1 ID %d parent 0%s", I=%d, ZONETYPE=ORDERED, DATAPACKING=POINT\\n'
>>> _ = open(t, "w").write(head + zl % (1, "", 2, ) + 2 * (row % 1) + zl % (1, " image 1", 1) + row % 1
...                        + zl % (2, "", 1) + row % 2)
>>> [(z["id"], z["image"], z["I"], len(z["rows"])) for z in parcel_zones(t)]
[(1, 0, 2, 2), (1, 1, 1, 1), (2, 0, 1, 1)]
>>> _ = open(t, "w").write(head + zl % (1, "", 2) + row % 1)
>>> parcel_zones(t)
Traceback (most recent call last):
    ...
ValueError: trajectories-A.dat: zone 1 (Mat A Group  1 ID 1 parent 0) holds 1 rows, its I= says 2
>>> _ = open(t, "w").write(head + zl % (2, "", 1) + row % 2 + zl % (1, "", 1) + row % 1)
>>> parcel_zones(t)
Traceback (most recent call last):
    ...
ValueError: trajectories-A.dat: zone 2 (Mat A Group  1 ID 1 parent 0) after group 1 ID 2 image 0: not in parcel order
"""
import os
import re

INT = re.compile(r"[+-]?\d+$")
PARCEL = re.compile(r"Mat (\S+) Group\s*(\d+) ID (\d+) parent (\d+)(?: image (\d+))?$")


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
                i = re.search(r"\bI\s*=\s*(\d+)", s[m.end():] if m else s)
                cur = {"title": m.group(1) if m else "", "names": list(names or []), "rows": [],
                       "I": int(i.group(1)) if i else None}
                out.append(cur)
                continue
            if low.startswith("title"):
                continue
            if names is None:
                raise ValueError("%s:%d: a data row before the variables= header" % (os.path.basename(path), n))
            if cur is None:
                cur = {"title": "", "names": list(names), "rows": [], "I": None}
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


def group(title):
    """The group number of a zone title, a group's or a parcel's; None without one."""
    m = re.search(r"Group\s*(\d+)", title)
    return int(m.group(1)) if m else None


def parcel(title):
    """(material, group, ID, parent, image) of a trajectory zone's title, image 0 for a parcel's first image and
    parent 0 for an injected parcel; None for another title."""
    m = PARCEL.match(title)
    if m is None:
        return None
    mat, g, i, p, k = m.groups()
    return mat, int(g), int(i), int(p), int(k or 0)


def parcel_zones(path):
    """The zones of a trajectory file, each with its title's fields (mat, group, id, parent, image). ValueError unless
    every zone is a parcel image's -- its title parses, I= equals its rows, every row carries the title's ID -- and the
    zones are in parcel order: by group, then parcel ID, a parcel's images together in ascending order."""
    out, prev = zones(path), None
    for n, z in enumerate(out, 1):
        where = "%s: zone %d (%s)" % (os.path.basename(path), n, z["title"])
        f = parcel(z["title"])
        if f is None:
            raise ValueError("%s: not a parcel image's title" % where)
        z.update(zip(("mat", "group", "id", "parent", "image"), f))
        if z["I"] != len(z["rows"]):
            raise ValueError("%s holds %d rows, its I= says %s" % (where, len(z["rows"]), z["I"]))
        bad = [r.get("ID") for r in z["rows"] if r.get("ID") != z["id"]]
        if bad:
            raise ValueError("%s: %d rows of another ID (%s)" % (where, len(bad), bad[0]))
        key = (z["group"], z["id"], z["image"])
        if prev is not None and key <= prev:
            raise ValueError("%s after group %d ID %d image %d: not in parcel order" % ((where,) + prev))
        prev = key
    return out


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
