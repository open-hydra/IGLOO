# wall-approach — the ord2 gas ghost ring must know its boundary conditions

Gates that a parcel drifting toward a solid wall is stopped by the gas instead of pushed
through it.

## What it covers

Under `gas-order = 2` the gas is sampled on a dual mesh whose boundary row straddles the
domain face: half of the first dual cell lies outside, on a ring of ghost nodes.
[fillGhostGradient](../../../src/lib/obj_block.f90) refills that ring at the start of every
`solve`, before any parcel reads it (`import_gas` rewrites only the interior): a linear
extrapolation, `2q1 − q2`, on every face, then a bc-aware pass on the gas-solid planes
(300/301); `fillGhostPartners` then fills the ghosts of connected faces.

The linear fill alone is wrong at a gas-solid plane in the one way a parcel can feel: the
sampled normal velocity in the boundary row becomes

    v_n = (v_ghost + v_1)/2 = (3 v_1 − v_2)/2 ,   not 0

so the gas keeps pushing the parcel into a plane it cannot itself cross. A bc-aware fill
written only at setup is inert: the per-sweep fill overwrites every ghost index before any
parcel reads it (measured: with the whole ord2 ghost ring poisoned to −9.9e30 at the end of
setup, 12 ord2 fixtures stayed byte-identical across 54 output files).

The bc-aware pass mirrors the **interior partner's** normal component into the ghost,

    v_g ← v_g − ((v_g + v_i)·n) n      ⇒   v_g·n = −v_i·n

keeping the tangential part from the linear fill. Mirroring the linear ghost instead
(`v_lin − 2(v_lin·n)n`) zeroes the plane value only for a uniform normal profile — see G2 in
[test_ghost_bc](../ghost_bc/test_ghost_bc.f90).

**Face 3 is tagged 301, not 300.** ATLAS writes 300 for the gas and 301 for every dispersed
phase of a `type = wall` patch, and IGLOO reads `<name>-bc.txt` — so a gas-solid wall reaches
IGLOO as 301. A mirror set of `{300}` alone would be inert on every generated case.

## Fixture

100 × 20 × 1 cells, 0.5 × 0.1 × 0.005 m (Δx = Δy = 5 mm), uniform gas `ρ 1, u 5, v −0.5,
T 350, γ 1.4, R 287, μ 2e-5, k 0.03`. Three 4 µm parcels are DB-injected at `x = 10 mm`,
`y = 36/41/46 mm`, `z = 2.5 mm`, at the gas velocity. `ρ_p = 2950` comes from
[common/properties.dat](../../common/properties.dat) column 3 — `[GPB-Phase1] rho` is an ATLAS
key IGLOO never reads.

`make_fixture.py` (committed, run once) writes both `INPUT/solfile.tec` and `INPUT/bc.txt`:
face 3 a 301 wall, faces 5/6 300, every other face 400 (no inlet: the parcels are DB-injected).
`tools/make_box_case.py` cannot tag face 3 as a 301 wall (its face 3 is 400, 300 or 201), and
planar-slab's `bc.txt` is ATLAS's, retagged by `retag_bc.py`. The injection heights sit off
every dual node plane (2.5 + 5k mm) and leave exit heights ≥ 4e-6 m (at 21/31/41 mm parcel 1
would exit at ~9e-9 m, 90× the ODE `atol` and below the `F12.6` print floor for some twenty
rows).

## Oracle

Above the wall row the gas is uniform, so each parcel drifts at `dy/dx = −0.1` and enters the
row (`|y| ≤ Δy/2`) at `x_s = 0.010 + 10 (y_0 − 0.0025)` = 0.345, 0.395, 0.445 m. Inside it the
sample is linear between the ghost `+V` and the first interior node `−V`, i.e.
`v_y(y) = −(2V/Δy) y` with `k = 200 s⁻¹`; with `τ_p = ρ_p d²/(18µ) = 1.31e-4 s` against
`1/(4k) = 1.25e-3 s` the parcel is overdamped by 9.5× and decays as `exp(−(k/U) x)`,
`k/U = 40 m⁻¹`, without ever crossing `y = 0`.

| check | asserts |
|---|---|
| W0 | ≥ 10 trajectory rows per parcel in the wall row — **non-vacuity**: without it W1/W2/W4 all pass on a run whose parcels never enter the row (a solfile written with `v = 0`). Measured 32/32/23 |
| W1 | no outloc row short of the outlet — i.e. zero wall hits |
| W2 | all three parcels exit at `x = L` |
| W3 | `y` non-increasing and ≥ 0 along every trajectory (ties allowed: the print floor) |
| W4 | the wall-row decay matches the tracer model to 5 % + 0.01. Measured residual **0.54 of the bound**; the excess is the finite-τ lag (`k τ_p = 0.026`), which makes the parcel fall slightly *faster* than the gas, not slower |
| W5 | the run took the planar single-layer path |

## Verified

**RED** on a solver whose ghost ring is the linear fill alone: all three parcels reach the
wall and are `gone` on contact, at `x = 0.370000 / 0.420000 / 0.470000` — the a-priori prediction
`x = 0.010 + 10 y_0` to six digits. W1 fires 3/3 and W4 at 25.5× the bound.
**GREEN** here: 0 hits, three outlet rows at `x = 0.5`, exit heights
`4.0e-6 / 3.3e-5 / 2.6e-4 m`, decay residual 0.54/0.54/0.51 of the bound.

The companion unit gate is [ghost_bc](../ghost_bc/test_ghost_bc.f90), which pins the faces,
edges, corners, the positivity guard and the ghost nodes of connected faces (101/201/102), the
last also end to end in `two-fam-bc-2blk-stretch-ord2` and `chimera-2blk-shift-ord2`.
