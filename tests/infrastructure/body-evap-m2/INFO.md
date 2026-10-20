# body-evap-m2 — the body-force accumulators of an evaporating parcel (ODE model 2)

**Purpose.** A body acceleration gives a parcel momentum and energy the gas never gave it, so the coupling source must
take that share out: each segment's deposit carries P += ṅ J g and E += ṅ W, with J = ∫m dt and W = ∫m g·v dt
integrated as extra ODE states for the models whose mass changes (2, 4 and 5; models 1 and 3 use a closed form). Model
2 keeps J and W at the last two slots, or W alone with the Eulerian output on, J then being the Eulerian mass moment.
`coupled-body` covers model 1 only.

**Fixture.** `INPUT/` is `evaporation/d2law-line`'s (symlink): the shared box, 25 drops of 30 µm at the inlet cell
centres, κ_v = 1, T₀ = 300 K, ρ_p = 8000 and a 100× c_p (125000 J/(kg K)) that freezes B_T, d2-law evaporation
(L_v = 2·10⁵). `input.ini` is d2law-line's with `body-accel = 0 −200 0`, `out-file = ALL`, `out-time = on` and
`time-end = 6.875e-3`: every parcel stops mid-cell 28, before the first one reaches the wall (7.5 ms), so its last row
is its final state and its exit record, code 0, is that state. `check.py` also runs the case with `out-file = S` in `ref-S/`, the
other slot layout.

**Oracle** (`check.py` → [`../../tools/body_source_check.py`](../../tools/body_source_check.py)). Momentum and energy
telescope along each path; the trajectory rows, which the source path does not write, give them:
- B1 ΣF_y = Σ_p ṅ[(m v_y)_first − (m v_y)_last] + g_y Σ_p ṅ∫m dt
- B2 ΣE = Σ_p ṅ[(m(h + |v|²/2))_first − (…)_last] + g_y Σ_p ṅ∫m v_y dt, h = c_p T

with ṅ = κ_ρ/(1 − κ_ρ) ρ_g U Δy Δz / m₀ (the injected number rate), the integrals by the trapezoid over the printed time
with a Richardson estimate. The enthalpy flux m c_p T is about 10⁸ times the body work here and the E13.6 print of m
alone is about 10⁻⁶ of it, so B2 values the mass the parcels lost by the mass source:
Σ_p ṅ m_last h_last = h_last (Σ_p ṅ m₀ − Σẇ), the 25 parcels sharing one final state within the integrator's spread
(h_last at their mean T, the spread in B2's tolerance) and B0 holding Σẇ to the rows.

**Gates** (each on both runs).
- W the witnesses: the body acceleration and the output fields in the log, the ODE size of the layout (15 with the
  Eulerian output, 10 without), one exit record per parcel, code 0 at its snapshot row's position and time, the
  injection state on every first row, every last row at t = time-end
  and in the snapshot; the parcels' final T within 5e-5 K of their mean and m within 1e-5 (B2's premise; measured 0).
- B0 Σẇ = Σ_p ṅ(m₀ − m_last) and ΣF_x = u_g Σ_p ṅ(m₀ − m_last) within the print of m; ΣF_z = 0.
- B1 within 1e-5 of the gross plus three Richardson estimates.
- B2 within the prints of T and of the kinetic term, the final states' spread and three Richardson estimates.
- B3 d² on the d2-law kernel 8 k_g/(c_p,g ρ_p) ln(1 + c_p,g(T_g − T)/L_v) along the measured T (the stagnant-film law
  does not see the slip), within the prints and the Richardson estimate.
- B4 non-vacuity: g_y Σ_p ṅ∫(m₀ − m) dt, what a closed form J = m₀ t would add, exceeds 100 B1 tolerances.
- B5 a regression band, no oracle: the two runs' `source.tec` agree cell by cell within 2.4e-6 of each variable's
  scale (10× the difference measured between the two layouts, which integrate one ODE under different error norms).

**Measured.** B0 0.41 of the m print; B1 6.5e-8 N (0.07 of its tolerance); B2 4.4e-4 W (0.60 — the T print, shared by
the 25 identical final states), the body work 1.08e-2 W being 15 tolerances; B3 0.47; B4 3167 B1 tolerances; B5 2.33e-7.

**RED.**
- Models 2 and 5 sent to the closed form of models 1 and 3 (`case(2,5)` of the deposit correction): B1 fails by
  3167 tolerances and B2 by 3.7 in both runs.
- The sign of W flipped in model 2's right-hand side: B2 fails by 28.8 tolerances in both runs.
- Half the d2-law rate: B3 fails on 700 rows (1.4e4 tolerances); B0–B2 hold (they close for any mass law).
