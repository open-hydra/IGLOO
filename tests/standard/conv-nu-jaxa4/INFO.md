# conv-nu-jaxa4 — the JAXA4 heat law at Mach > 0

**Purpose.** `conv-nu`'s constant-slip state under `heat = JAXA4`,
Nu = [(2 + 0.654 Re^½ Pr^⅓)^−1 + 3.42 Ma/(Re Pr)]^−1 as `[Shim06]` eq. 50 prints it: the continuum base and the
rarefaction term both enter the temperature history, through the solver's Re, Pr and slip Mach number. NASA-SP-8039
Table II prints the same law with Pr^½ in its last row; the code follows Shimada (Pr^⅓, the exponent of every other
law of the table), and so does this oracle — with Pr^½ the case would move by −6.6e-3 of Nu, up to 0.73 K.

**Fixture.** `INPUT/` is `conv-nu`'s (symlink): the shared box, `body-accel = 0 −200 0`, 25 parcels injected at their
terminal velocity (κ_v = 1.00033, α = −0.02574 rad), κ_t = 0.5 (T₀ = 300 K in a 600 K gas), d = 11.89 µm, Stokes drag.
`input.ini` is `conv-nu`'s with `heat = JAXA4` and `out-time = on`.

**Oracle** (`check.py`, the shared [`../conv-nu/nu_law_check.py`](../conv-nu/nu_law_check.py)). Re = 0.204062,
Pr = 0.695423 and Ma = |v_∞|/√(γRT_g) = 5.2431e-4 are constant along every path (Ma/(Re Pr) = 3.6947e-3), so
T(t) = T_g + (T₀ − T_g) e^{−t/τ_T}, τ_T = c_p ρ_p d²/(6 Nu k_g), with Nu = 2.198900 from the printed form (2.261743 at
Ma = 0: the rarefaction term is 2.78 % of Nu) and t the printed parcel time.

**Gates.** G0 the injection state on every first row (t = 0, T = 300 K); G1 `conv-nu`'s guards (u = u_g, v = v_∞,
w = 0, d = d_p) on every row; G2 T on every row within 5e-7 + |dT/dt|·5e-9·t + 1e-8 K (the T print, the time print,
the integrator), 25 parcels verified; N1 the Ranz-Marshall curve and the law's own Ma = 0 curve each miss the rows by
more than 1e3 tolerances (the case sees the law and its Mach term).

**Measured.** G2 worst 0.961 of the tolerance, the print itself: no printed T moves between the ODE tolerances 1e-11
and 1e-13, and the same run at 1e-7 fails G2 at 1.86 of the tolerance. N1 1.9e6 (Ranz-Marshall), 2.9e6 (Ma = 0).

**RED.**
- The Mach number withheld from `heat()` in `interphase` (`heat(Re,Pr,0,…)`): G2 fails by 3.110 K (2.9e6
  tolerances), the rows then on the Ma = 0 curve within 0.93 of the tolerance; `conv-nu-jaxa2/-jaxa3` stay green.
- The base 0.654 → 0.645 in `heat_JAXA_4`: G2 fails by 0.171 K (1.6e5 tolerances), `test_heat_nu`'s HN6 by 1.3e-2.
