#!/usr/bin/env python3
"""INPUT/properties.dat for khrt-varcp: khrt-e2e's water (density 1000) with a varying cp on 1..1000 K,
cp = 4182 + 10 (T - 270) J/(kg K), and the enthalpy h(1 K) = cp(1) then its trapezoid sum (exact for a linear cp:
h = cp(1) + 4182 (T - 1) + 5 ((T - 270)^2 - 269^2)). Run once; the output is committed."""
T0, T1, RHO = 1, 1000, 1000.0


def cp(T):
    return 4182.0 + 10.0 * (T - 270)


def h(T):
    return cp(T0) + 4182.0 * (T - T0) + 5.0 * ((T - 270) ** 2 - (T0 - 270) ** 2)


with open("INPUT/properties.dat", "w") as g:
    g.write('TITLE = "Mass Thermodynamic Properties"\n')
    g.write('VARIABLES = "Temperature", "Cp", "Density", "Enthalpy"\n')
    g.write('ZONE T="A"\n')
    g.write("I=%d, F=POINT\n" % (T1 - T0 + 1))
    for T in range(T0, T1 + 1):
        g.write("%.1f %.6f %.6f %.6f\n" % (T, cp(T), RHO, h(T)))
print("INPUT/properties.dat written: %d..%d K" % (T0, T1))
