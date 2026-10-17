#!/usr/bin/env python3
"""INPUT/properties.dat for d2law-varcp: d2law's material (density 2950) with a varying cp on 250..900 K,
cp = 1000 + 2 (T - 250) J/(kg K), and the enthalpy its trapezoid sum from h(250 K) = 2.5e5 J/kg (exact for a linear
cp: h = 1000 T + (T - 250)^2), the relative datum. Run once; the output is committed."""
T0, T1, RHO = 250, 900, 2950.0
with open("INPUT/properties.dat", "w") as g:
    g.write('TITLE = "Mass Thermodynamic Properties"\n')
    g.write('VARIABLES = "Temperature", "Cp", "Density", "Enthalpy"\n')
    g.write('ZONE T="A"\n')
    g.write("I=%d, F=POINT\n" % (T1 - T0 + 1))
    for T in range(T0, T1 + 1):
        g.write("%.1f %.6f %.6f %.6f\n" % (T, 1000.0 + 2.0 * (T - T0), RHO, 1000.0 * T + (T - T0) ** 2))
print("INPUT/properties.dat written: %d..%d K" % (T0, T1))
