#!/usr/bin/env python3
"""INPUT/properties.dat for evap-breakup-box-varcp: evap-breakup-box's material (density 2950) with a varying cp on
1..1000 K, cp = 1250 + (T - 300) J/(kg K), and the enthalpy h(1 K) = cp(1) then its trapezoid sum (exact for a linear
cp: h = cp(1) + 1250 (T - 1) + ((T - 300)^2 - 299^2)/2). Run once; the output is committed."""
T0, T1, RHO = 1, 1000, 2950.0


def cp(T):
    return 1250.0 + (T - 300)


def h(T):
    return cp(T0) + 1250.0 * (T - T0) + 0.5 * ((T - 300) ** 2 - (T0 - 300) ** 2)


with open("INPUT/properties.dat", "w") as g:
    g.write('TITLE = "Mass Thermodynamic Properties"\n')
    g.write('VARIABLES = "Temperature", "Cp", "Density", "Enthalpy"\n')
    g.write('ZONE T="A"\n')
    g.write("I=%d, F=POINT\n" % (T1 - T0 + 1))
    for T in range(T0, T1 + 1):
        g.write("%.1f %.6f %.6f %.6f\n" % (T, cp(T), RHO, h(T)))
print("INPUT/properties.dat written: %d..%d K" % (T0, T1))
