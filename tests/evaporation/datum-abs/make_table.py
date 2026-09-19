#!/usr/bin/env python3
"""INPUT/properties.dat for datum-abs: tc-box's constant-cp table with the Enthalpy column shifted
to an ABSOLUTE datum, h(T) = cp*T + H_OFF, and tagged Enthalpy_abs (what ATLAS GPB writes for a
fixed-cp material given [GPB-Phase*] h0 = cp*298.15 + H_OFF). Run once; the output is committed."""
H_OFF = -1.7113460e7   # J/kg: the water value of hydra's evap-box cases (h0 = -15866000, cp = 4184)
src = "../tc-box/INPUT/properties.dat"
with open(src) as f, open("INPUT/properties.dat", "w") as g:
    for line in f:
        t = line.split()
        if line.startswith("VARIABLES"):
            g.write(line.replace('"Enthalpy"', '"Enthalpy_abs"'))
        elif len(t) == 4 and t[0][0].isdigit():
            g.write(f"{t[0]} {t[1]} {t[2]} {float(t[3]) + H_OFF:.6f}\n")
        else:
            g.write(line)
print("INPUT/properties.dat written with H_OFF =", H_OFF)
