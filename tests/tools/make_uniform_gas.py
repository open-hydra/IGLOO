#!/usr/bin/env python3
"""Constantize a known-good IGLOO gas Tecplot file into a spatially-UNIFORM field.

Reads a real `solfile.tec` (BLOCK packing, 1 value/line), keeps the header and the
nodal mesh coordinates (X,Y,Z) byte-for-byte, and overwrites every cell-centered
variable block with a single constant chosen by VARIABLE NAME (IGLOO's gas import
in allocation.f90 maps by name, not column). Particle-side variables (rho_p, u_p,
...) are zeroed -- IGLOO never reads them.

The output is therefore guaranteed format/mesh/BC-consistent with the source
(same reader, same mesh), differing only in that the gas state is uniform. That
makes the per-particle trajectory an EXACT Stokes/relaxation problem with a
closed-form solution, independent of the (real, curved) mesh geometry.

Usage:
    make_uniform_gas.py <src solfile.tec> <dst solfile.tec>

Constants are defined in GAS below (SI units). Edit there to retune a case.
"""
import re
import sys

# Uniform gas state (SI). Names match the substring/exact rules in
# src/lib/allocation.f90::import_gas.
GAS = {
    "Roi": 1.2,      # density        [kg/m^3]   (matches 'Roi(' or 'rho(')
    "U":   10.0,     # x-velocity     [m/s]      (exact 'U')
    "V":   0.0,      # y-velocity     [m/s]      (exact 'V')
    "W":   0.0,      # z-velocity     [m/s]      (exact 'W')
    "P":   101325.0, # pressure       [Pa]       (NOT read by IGLOO)
    "T":   600.0,    # temperature    [K]        (exact 'T')
    "MIL": 1.8e-5,   # laminar visc.  [Pa s]     (matches 'MIL')
    "KL":  0.026,    # conductivity   [W/m/K]    (matches 'KL')
    "GAM": 1.4,      # gamma          [-]        (matches 'GAM')
    "R":   287.0,    # gas constant   [J/kg/K]   (matches 'R', post-density)
    "MIT": 1.8e-5,   # total visc.    [Pa s]     (matches 'MIT'); = MIL so drag mu unambiguous
}


def const_for(varname):
    """Return the constant for a cell-centered variable, by the same matching
    logic IGLOO uses; particle vars (and anything unmatched) -> 0."""
    v = varname.strip().strip('"')
    if "Roi(" in v or "rho(" in v:
        return GAS["Roi"]
    if v in ("U", "u"):
        return GAS["U"]
    if v in ("V", "v"):
        return GAS["V"]
    if v in ("W", "w"):
        return GAS["W"]
    if v == "T":
        return GAS["T"]
    if "MIL" in v or "mil" in v:
        return GAS["MIL"]
    if "MIT" in v or "mit" in v:
        return GAS["MIT"]
    if "KL" in v or "kl" in v:
        return GAS["KL"]
    if "GAM" in v:
        return GAS["GAM"]
    if v in ("R",):            # the exact gas-constant column
        return GAS["R"]
    if v in ("P", "p"):
        return GAS["P"]
    return 0.0                  # particle-side / unread variables


def main():
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    src, dst = sys.argv[1], sys.argv[2]

    with open(src) as f:
        lines = f.readlines()

    # --- locate header pieces dynamically (header length varies: DATAPACKING and
    #     VARLOCATION may share a line or be split across two) ---
    var_idx = next(i for i, l in enumerate(lines) if "VARIABLES" in l.upper())
    zone_idx = next(i for i, l in enumerate(lines)
                    if re.search(r"\bI\s*=\s*\d+", l) and re.search(r"\bK\s*=\s*\d+", l))
    loc_idx = next((i for i, l in enumerate(lines) if "VARLOCATION" in l.upper()), None)
    pack_idx = next((i for i, l in enumerate(lines) if "DATAPACKING" in l.upper()), zone_idx)
    header_len = (loc_idx if loc_idx is not None else pack_idx) + 1

    varnames = re.findall(r'"([^"]*)"', lines[var_idx])
    nvar = len(varnames)

    m = re.search(r"I\s*=\s*(\d+).*?J\s*=\s*(\d+).*?K\s*=\s*(\d+)", lines[zone_idx])
    I, J, K = int(m.group(1)), int(m.group(2)), int(m.group(3))

    loc_text = lines[loc_idx] if loc_idx is not None else ""
    mloc = re.search(r"\[\s*(\d+)\s*-\s*(\d+)\s*\]\s*=\s*CELLCENTERED", loc_text, re.I)
    if mloc:
        cc_lo, cc_hi = int(mloc.group(1)), int(mloc.group(2))
    else:                               # no VARLOCATION => all vars nodal
        cc_lo, cc_hi = nvar + 1, nvar

    nnode = I * J * K
    ncell = max(I - 1, 1) * max(J - 1, 1) * max(K - 1, 1)

    # value layout: vars 1..(cc_lo-1) are nodal (nnode each), cc_lo..cc_hi are
    # cell-centered (ncell each). Any trailing nodal vars would follow, but here
    # all of 4..nvar are cell-centered.
    counts = []
    for vi in range(1, nvar + 1):
        counts.append(ncell if cc_lo <= vi <= cc_hi else nnode)
    expected = header_len + sum(counts)
    if expected != len(lines):
        sys.exit(f"[ERROR] line-count mismatch: header+values={expected} "
                 f"but file has {len(lines)} lines (header_len={header_len} "
                 f"I={I} J={J} K={K} nnode={nnode} ncell={ncell} nvar={nvar})")

    out = lines[:header_len]
    pos = header_len
    for vi in range(1, nvar + 1):
        n = counts[vi - 1]
        block = lines[pos:pos + n]
        if cc_lo <= vi <= cc_hi:
            c = const_for(varnames[vi - 1])
            const_str = f"  {c:.15E}\n"
            out.extend([const_str] * n)        # overwrite cell var -> constant
        else:
            out.extend(block)                  # keep nodal mesh coords verbatim
        pos += n

    with open(dst, "w") as f:
        f.writelines(out)

    print(f"[ok] {dst}: I={I} J={J} K={K} nnode={nnode} ncell={ncell} "
          f"nvar={nvar} cellcentered=[{cc_lo}-{cc_hi}]")
    cc = [(varnames[vi - 1].strip(), const_for(varnames[vi - 1]))
          for vi in range(cc_lo, cc_hi + 1)]
    print("      cell-centered constants:")
    for name, val in cc:
        print(f"        {name:8s} = {val}")


if __name__ == "__main__":
    main()
