# sym-plane-hop — parcels bouncing in place on a symmetry plane: the reflection's crossing

**Purpose.** A reflection off a symmetry plane (code 300) takes its crossing — position, velocity and time — from the
step that crossed the face: the cubic through the step's two end states, its first crossing of the face plane
(`obj_bc.f90::hopCrossing`). A parcel pressed onto the plane rebounds in hops that grow shorter than the solver's steps:
a step then holds a whole hop, leaving the plane and coming back below it, with a chord shorter than the refinement's
`eps` = 1e-8 m. The crossing taken on the straight chord there has the step's start velocity, pointing away from the
plane, so its mirror image sends the parcel into the plane, and the time rolled back from the step's end over the
chord's mean speed can fall before the step's start. This case bounces three parcels in place against a closed form.

**Fixture.** `stuck-cell`'s `solfile.tec` (the box `tools/make_box_case.py --sym-face 5 --no-inlet --ug 0.0` writes:
drag-stokes's 60×5×5 cells over 0.15 × 0.05 × 0.05 m, the gas at rest, T = 600 K, μ = 1.8e-5 Pa s) and
`sym-plane-skate`'s `bc.txt` (face 5, z = 0, a symmetry plane; no inlet). Three `[IGLOO-BC]` streams at x = 0.004 m,
y = 0.015, 0.025, 0.035 m, z₀ = 5e-5 m, at rest along the plane and falling at w₀ = −0.02, −0.05, −0.08 m/s;
d = 11.89 µm, ρ_p = 2950 kg/m³ (τ = ρ_p d²/(18 μ) = 1.287184e-3 s), `drag = Stokes`, `heat = NoHeat` at the gas
temperature (no exchange either way), `body-accel = 0 0 -9.81`, SDIRK4 at its default tolerances, gas order 1,
`print-dtime = 1e-300` (a row at every segment end), `out-time = on`. Every segment after the first ends in a reflection
in the injection cell (10 mm high, hops below 0.2 mm), and the stuck-bounce guard ends each parcel at its 100th
reflection, within 5e-5 m of its injection point (under the guard's 1e-4 m).

**Oracle** (in `check.py`). Between impacts the normal motion is closed-form under linear drag and gravity,
w(t) = −gτ + (w₀ + gτ) e^(−t/τ), z(t) = z₀ − gτ t + τ (w₀ + gτ)(1 − e^(−t/τ)); each impact is a root of z = 0
(bisection to machine precision) and each rebound the impact speed reflected, 1e-12 m above the plane as the
reflection places it. The 100th reflection at t₁₀₀ = 1.9750386654e-2, 2.0515580573e-2 and 2.3135439956e-2 s; the
hops are shorter than 8e-5 s from the 49th–50th reflection on, 3.9e-5 s at the 100th.

**Gates** (`check.py`).
- H0 (precondition): three exit records, code 0, on the plane (z ≤ 1e-5 m) at their injection (x, y); one
  `stuck in cell` line per stream; no `no net progress`, `hit outer maxIter`, `Inner loop >`, `Run_ODESolver err=` or
  `non-finite` line; the run ends normally.
- H1 the reflection times: |t_exit − t₁₀₀|/t₁₀₀ ≤ 5e-3 per stream. The crossing on the step's cubic is exact under
  constant acceleration and third-order otherwise: over hops 40–100 and any step of 1–2 hops its error is at most
  4.53e-4 of a hop in time and 8.92e-4 in the rebound speed, which accumulated to the 100th reflection gives at most
  3.31e-3 even with every hop from the 40th a single step of three hops (an emulation of the rule on whole-hop steps).
  Measured +4.35e-6, +3.83e-6, +3.66e-6.
- H2 the time never falls back: each trajectory has 1 + 100 rows (the injection and one per segment end, each a
  reflection, the 100th the last) with strictly increasing t. Measured 101 rows each; every row's time within 2.2e-3
  hops of the closed form's reflection (printed, not gated).
- H3 the rebound leaves the plane: w > 0 on every row after the injection row. Measured on all 100.

**RED.** The crossing on the straight chord (the rule kept for a step without usable data): H1 −8.79e-2, −7.64e-2,
−7.67e-2; H2 82, 92 and 82 rows (a segment whose time falls before the last row writes none); H3 14, 11 and 12 rows
with w < 0, the first at rows 64, 60 and 60. Into-the-plane rebounds (the crossing velocity pointing back into the
domain) number 17, 11 and 16; the rows the time skips also follow outward crossings whose rollback exceeds the step
(up to 9.3e-5 s, more than two hops).
