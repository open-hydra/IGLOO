# body-evap-m4 — the body-force accumulators of an evaporating parcel under ODE model 4

**Purpose.** As [`body-evap-m2`](../body-evap-m2/INFO.md), for model 4 (evaporation with an ODE breakup rate), whose
states carry the stream's number rate: J = ∫ṅ m dt and W = ∫ṅ m g·v dt. With the Eulerian output on, J is the
Eulerian mass moment at neq − 2 (neq − 1 holds the number rate); with it off, J and W sit at the last two slots. The
deposit correction reads J from the slot of the layout in force.

**Fixture.** `INPUT/` is `evaporation/d2law-brk-dormant`'s (symlink, itself `d2law`'s): the shared box, 25 drops of
30 µm, κ_v = 1, T₀ = 300 K, ρ_p = 2950, c_p = 1250, d2-law evaporation (L_v = 2·10⁵) with the Reitz-Diwakar breakup rate
switched on (σ = 0.072, μ_p = 10⁻³). `input.ini` is d2law-brk-dormant's with `body-accel = 0 −200 0`, `out-file = ALL`,
`out-time = on` and `time-end = 6.875e-3` (every parcel stops mid-cell 28, before any reaches a wall). `check.py` also
runs the case with `out-file = S` in `ref-S/`.

**Oracle and gates** (`check.py` → [`../../tools/body_source_check.py`](../../tools/body_source_check.py)): W, B0–B5 as
in `body-evap-m2` (ODE size 17 with the Eulerian output, 11 without; B5's band 1.9e-6), plus
- D the breakup stays dormant under the gravity slip: We < 6 and We/Re^½ < 0.5 on every row, and the Eulerian number
  total Σ n_p V equals 25 ṅ time-end (the number rate never moved).

**Measured.** B0 0.06 of the m print; B1 7.9e-7 N (0.11 of its tolerance); B2 5.3e-6 W (0.14), the body work 6.39e-3 W
being 166 tolerances; B3 0.43; B4 913 B1 tolerances; D We ≤ 2.4e-4, We/Re^½ ≤ 2.2e-4, Σ n_p V exact to 13 digits;
B5 1.82e-7.

**RED.**
- The deposit correction reading J at neq − 1 whatever the layout (the Eulerian case of `jSlot` dropped): the
  out-file = ALL run takes the number-rate integral for J and B1 fails by 5.1e8 N; `ref-S` stays green; B5 fails
  (1.0 of the scale).
- The sign of W flipped in model 4's right-hand side: B2 fails by 332 tolerances in both runs.
- Half the d2-law rate: B3 fails on 700 rows (1.2e4 tolerances).
