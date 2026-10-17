# conv-nu-jaxa3 — the JAXA3 heat law at constant slip

**Purpose.** `conv-nu`'s constant-slip state under `heat = JAXA3`, Nu = 2 + 0.459 Re^0.55 Pr^⅓
(`[SP8039]` Table II, ref. 42): the law enters the temperature history end to end, through the solver's Re and Pr.

**Fixture.** `INPUT/` is `conv-nu`'s (symlink): the shared box, `body-accel = 0 −200 0`, 25 parcels injected at their
terminal velocity (κ_v = 1.00033, α = −0.02574 rad), κ_t = 0.5 (T₀ = 300 K in a 600 K gas), d = 11.89 µm, Stokes drag.
`input.ini` is `conv-nu`'s with `heat = JAXA3` and `out-time = on`.

**Oracle** (`check.py`, the shared [`../conv-nu/nu_law_check.py`](../conv-nu/nu_law_check.py)). Re = 0.204062,
Pr = 0.695423 and the slip Mach number Ma = 5.2431e-4 are constant along every path, so
T(t) = T_g + (T₀ − T_g) e^{−t/τ_T}, τ_T = c_p ρ_p d²/(6 Nu k_g), with Nu = 2.169667 from the printed form (typed in the
check, not read from the code) and t the printed parcel time.

**Gates.** G0 the injection state on every first row (t = 0, T = 300 K); G1 `conv-nu`'s guards (u = u_g, v = v_∞,
w = 0, d = d_p) on every row; G2 T on every row within 5e-7 + |dT/dt|·5e-9·t + 1e-8 K (the T print, the time print,
the integrator), 25 parcels verified; N1 the Ranz-Marshall curve misses the rows by more than 1e3 tolerances.

**Measured.** G2 worst 0.972 of the tolerance, the print itself: no printed T moves between the ODE tolerances 1e-11
and 1e-13, and the same run at 1e-7 fails G2 at 1.94 of the tolerance. N1 3.3e6.

**RED.** The coefficient 0.459 → 0.495 in `heat_JAXA_3`: G2 fails by 0.675 K (6.4e5 tolerances) and
`test_heat_nu`'s HN5 by 7.5e-2. The Mach number withheld from `heat()` leaves the case green: JAXA3 is Mach-free.
