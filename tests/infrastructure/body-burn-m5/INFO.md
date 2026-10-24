# body-burn-m5 — the body-force accumulators of a burning parcel (ODE model 5)

**Purpose.** As [`body-evap-m2`](../body-evap-m2/INFO.md), for model 5 (the Beckstead d^n burn law), whose slots are
model 2's: J = ∫m dt and W = ∫m g·v dt at the last two slots, or W alone with the Eulerian output on.

**Fixture.** `INPUT/` is `combustion/burn-box`'s (symlink): the shared box, 25 aluminium drops of 30 µm, κ_v = 1,
T₀ = T_g = 600 K, ρ_p = 2950, c_p = 1250, `combustion=Beckstead` with K-burn 4.5e-7, n 1.8, X-eff 0.5, T-ign 550 K,
q-comb 10⁶ J/kg. `input.ini` is burn-box's with `body-accel = 0 −200 0`, `out-file = ALL`, `out-time = on` and
`time-end = 6.875e-3` (every parcel stops mid-cell 28, before any reaches a wall). `check.py` also runs the case with
`out-file = S` in `ref-S/`. burn-box's own trajectory oracle does not apply: the slip changes the drop's heat.

**Oracle and gates** (`check.py` → [`../../tools/body_source_check.py`](../../tools/body_source_check.py)): W, B0–B2,
B4 and B5 as in `body-evap-m2` (ODE size 15 with the Eulerian output, 10 without; B5's band 1.0e-6); no B3.

**Measured.** B0 0.65 of the m print; B1 1.1e-7 N (0.09 of its tolerance); B2 6.9e-7 W (0.02), the body work
8.53e-3 W being 240 tolerances; B4 2669 B1 tolerances; B5 9.26e-8.

**RED.**
- Models 2 and 5 sent to the closed form of models 1 and 3 (`case(2,5,8)` of the deposit correction): B1 fails by
  2669 tolerances and B2 by 66 in both runs.
- The sign of W flipped in model 5's right-hand side: B2 fails by 479 tolerances in both runs.
