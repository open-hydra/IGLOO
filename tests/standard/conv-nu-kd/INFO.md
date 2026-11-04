# conv-nu-kd — the Kavanau-Drake heat law at Mach > 0

**Purpose.** `conv-nu`'s constant-slip state under `heat = Kavanau-Drake`,
Nu = N₀/(1 + 3.42 (Ma/(Re Pr)) N₀), N₀ = 2 + 0.459 Re^0.55 Pr^0.33 (`[SP8039]` Table II, ref. 72; `[Shim06]` eq. 49,
both with the exponent 0.33 printed): the continuum base and the rarefaction term both enter the temperature history,
through the solver's Re, Pr and slip Mach number.

**Fixture.** `INPUT/` is `conv-nu`'s (symlink): the shared box, `body-accel = 0 −200 0`, 25 parcels injected at their
terminal velocity (κ_v = 1.00033, α = −0.02574 rad), κ_t = 0.5 (T₀ = 300 K in a 600 K gas), d = 11.89 µm, Stokes drag.
`input.ini` is `conv-nu`'s with `heat = Kavanau-Drake` and `out-time = on`.

**Oracle** (`check.py`, the shared [`../conv-nu/nu_law_check.py`](../conv-nu/nu_law_check.py)). Re = 0.204062,
Pr = 0.695423 and Ma = |v_∞|/√(γRT_g) = 5.2431e-4 are constant along every path (Ma/(Re Pr) = 3.6947e-3), so
T(t) = T_g + (T₀ − T_g) e^{−t/τ_T}, τ_T = c_p ρ_p d²/(6 Nu k_g), with Nu = 2.111966 from the printed form (2.169873 at
Ma = 0: the rarefaction term is 2.67 % of Nu) and t the printed parcel time.

**Gates.** G0 the injection state on every first row (t = 0, T = 300 K); G1 `conv-nu`'s guards (u = u_g, v = v_∞,
w = 0, d = d_p) on every row; G2 T on every row within 5e-7 + |dT/dt|·5e-9·t + 1e-8 K (the T print, the time print,
the integrator), 25 parcels verified; N1 the Ranz-Marshall curve and the law's own Ma = 0 curve each miss the rows by
more than 1e3 tolerances (the case sees the law and its Mach term).

**Measured.** G2 worst 0.939 of the tolerance, the print itself: no printed T moves between the ODE tolerances 1e-11
and 1e-13; the same run fails G2 at 1e-8 (1.04) and 1e-7 (1.87). N1 6.1e6 (Ranz-Marshall), 2.8e6 (Ma = 0).

**RED.**
- The Mach number withheld from `heat()` in `interphase`: G2 fails by 2.983 K (2.8e6 tolerances).
- The exponent 0.33 → 1/3 in `heat_Kavanau_Drake` (9.2e-5 of Nu): G2 fails by 1.02e-2 K (9.6e3 tolerances),
  `test_heat_nu`'s HN7 by 6.2e-3 — the case keeps the printed 0.33 from being "cleaned up" silently.
- The rarefaction coefficient 3.42 → 3.24: G2 fails by 0.155 K (1.5e5 tolerances), HN7 by 5.3e-2, HN8 by 2.9e13
  of its bound.
