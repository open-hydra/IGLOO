#!/usr/bin/env python3
"""The two-zone properties tables of the oxide-cap cases: zone A the burning metal (burn-box's Cp 1250 and Density
2950), zone OX its combustion product, liquid alumina (Cp 1887.6, Density 3950). Run once; the outputs are committed.

  python3 make_table.py                      INPUT/properties.dat: Enthalpy_abs, A = 1250 T - 372687.5 (zero at
                                             298.15 K), OX = 1887.6 T - 1.723894e7 (JANAF's liquid line at 2327 K)
  python3 make_table.py varcp <path>         the same, OX's Cp rising 0.01 J/kg/K per K about 2327 K (its enthalpy
                                             the integral of it), ten rows (300-309 K): a refusal fixture
  python3 make_table.py relative <path>      Enthalpy = Cp T for both zones (the relative datum), ten rows
"""
import sys

CP_A, RHO_A, HOFF_A = 1250.0, 2950.0, -372687.5
CP_OX, RHO_OX, HOFF_OX = 1887.6, 3950.0, -1.723894e7


def zone(name, rows):
    out = ['ZONE T="%s"' % name, "I=%d, F=POINT" % len(rows)]
    out += ["%.1f %.6f %.6f %.6f" % r for r in rows]
    return out


def main():
    variant = sys.argv[1] if len(sys.argv) > 1 else "abs"
    path = sys.argv[2] if len(sys.argv) > 2 else "INPUT/properties.dat"
    T = [float(k) for k in (range(1, 5001) if variant == "abs" else range(300, 310))]
    if variant == "relative":
        a = [(t, CP_A, RHO_A, CP_A * t) for t in T]
        ox = [(t, CP_OX, RHO_OX, CP_OX * t) for t in T]
        head = '"Temperature", "Cp", "Density", "Enthalpy"'
    else:
        a = [(t, CP_A, RHO_A, CP_A * t + HOFF_A) for t in T]
        if variant == "varcp":
            ox = [(t, CP_OX + 0.01 * (t - 2327.0), RHO_OX, CP_OX * t + 0.005 * (t - 2327.0) ** 2 + HOFF_OX) for t in T]
        else:
            ox = [(t, CP_OX, RHO_OX, CP_OX * t + HOFF_OX) for t in T]
        head = '"Temperature", "Cp", "Density", "Enthalpy_abs"'
    lines = ['TITLE = "Mass Thermodynamic Properties"', "VARIABLES = " + head] + zone("A", a) + zone("OX", ox)
    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")
    print("%s written (%s, %d rows per zone)" % (path, variant, len(T)))


if __name__ == "__main__":
    main()
