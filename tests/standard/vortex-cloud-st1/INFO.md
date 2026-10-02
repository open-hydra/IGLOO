# vortex-cloud — a particle cloud in a frozen solid-body vortex, exact per parcel

**Purpose.** It is the suite's only case where the gas rotates in the plane of the motion, so the full two-dimensional coupling of the drag equation — every cross term between x and y — is gated, and the only one that sweeps the Stokes number: four directories, `vortex-cloud-st0p01`, `-st0p1`, `-st1` (this one, which holds `INPUT/` and `check.py`) and `-st10`, differing only in the particle diameter.

**Physics.** In the complex plane `z = x + i y` the gas is `u_g = i Ω z` (Ω = 1 rad/s). Stokes
drag with a uniform diameter gives, along a trajectory,

    τ z'' + z' = i Ω z,     λ± = (−1 ± sqrt(1 + 4 i Ω τ)) / (2τ),

and a parcel seeded at the local gas velocity, `z'(0) = i Ω z0`, moves as `z(t) = C(t) z0` with
`C = A e^{λ+ t} + B e^{λ− t}`, `A + B = 1`, `A λ+ + B λ− = i Ω`. One complex factor carries the whole
cloud: `arg C` lags `Ω t` and `|C|` exceeds 1 (inertia throws the parcels outwards). At a quarter
turn of the gas, `t = QUARTER = π/2`:

| St | d [m] | \|C\| | lag Ω t − arg C [rad] |
|---|---|---|---|
| 0.01 | 5.692e-5 | 1.015723 | 0.000310 |
| 0.1 | 1.800e-4 | 1.151551 | 0.025864 |
| 1 | 5.692e-4 | 1.617416 | 0.345426 |
| 10 | 1.800e-3 | 1.828830 | 0.538470 |

**Fixture** (`tests/tools/make_vortex_case.py` regenerate with
`python3 tests/tools/make_vortex_case.py tests/standard/vortex-cloud-st1 1 --input` and the other three without `--input`). Gas: the square `[−2, 2]²`, 96 × 96 cells (h = 1/24), one cell of thickness h across z (the planar single-layer path), cell-centred `u = −Ω y`, `v = Ω x`, ρ 1.2, T 300, μ 1.8e-5 (`mil` and `mit`), γ 1.4, R 287.05 — linear, so the second-order rebuild samples it exactly (`gas-order = 2`). Every face is 400. Material: `INPUT/properties.dat` at ρ_p 1000, cp 900. Cloud: a lattice of spacing h/2 centred on (0.6, 0) with half-spacing offsets (its weighted centroid is 0.6 to
4.4e-16): 1672 DB parcels, every seed at least 0.15 h from a node or dual line, seeded at `(−Ω y0, Ω x0)`, T 300, `mdot = ρ_p0(z0) (h/2)² h`. Stokes drag, `heat = NoHeat`, SDIRK4 at 1e-11, `time-end = QUARTER`, `out-time = on`, `out-file = S`. Each list of 1672 values is one 37.6 kB INI line; FiNeR's `count_values` reads 1672 on each.

**Gates** (`check.py`; the oracle is the closed form above, cross-checked by RK4):

| id | assertion |
|---|---|
| G0' | the input's seed lists are the generator's; all 1672 IDs start from their seed state at t = 0 (F12.6 half-ULP); `snapshot-A.dat` holds every parcel at t = QUARTER (1e-8 s); no exit record; the log reports the planar single-layer path and no lost, stuck or stalled parcel |
| G1 | `C(0) = 1`, `C'(0) = i Ω`; `Im(C'/C) > 0` and `d|C|/dt > 0` on [1e-3, QUARTER]; the closed form against RK4 of `τ z'' + z' = i Ω z` at QUARTER to 1e-12 |
| G2' | at every time-stamped row and snapshot record, `|z − C(t) z0| ≤ 3e-6` m and `|v − C'(t) z0| ≤ 3e-6` m/s |

Printed, not gated: the cloud's own `|C|`, lag and velocity slope from the snapshot (weights mdot
by ID) beside the closed form.

**Measured.** All four pass; the G2' residual is the F12.6 print floor (√2 · 5e-7 = 7.1e-7):

| St | rows | worst \|z − C z0\| [m] | worst \|v − C' z0\| [m/s] | cloud \|C\| / lag at QUARTER |
|---|---|---|---|---|
| 0.01 | 55873 | 7.01e-7 | 7.06e-7 | 1.015723 / 0.000310 |
| 0.1 | 58481 | 6.99e-7 | 7.05e-7 | 1.151551 / 0.025864 |
| 1 | 57008 | 6.98e-7 | 7.04e-7 | 1.617416 / 0.345426 |
| 10 | 54107 | 6.69e-7 | 7.06e-7 | 1.828830 / 0.538470 |

The tolerance sits 4.2× above the floor. Wall time ≈ 15 s per directory at OMP 5. Under MPI (4 ranks
× 2 threads) the St = 1 snapshot and trajectory rows are the serial ones as sorted multisets.

**RED runs** (deliberately wrong inputs; worst residual of G2' or the premise that refuses):

| id | wrong input | St 0.01 | St 0.1 | St 1 | St 10 |
|---|---|---|---|---|---|
| R1 | `drag = NoDrag` (straight flight, \|C\| = sqrt(1 + Ω²t²) = 1.862096 at QUARTER) | 1.23 m | 1.13 m | 0.488 m | 6.67e-2 m |
| R2 | particle density 2000 (τ doubled, input diameter unchanged) | 1.69e-2 m | 0.148 m | 0.201 m | 3.27e-2 m |
| R3 | seeded at rest (`up = vp = 0`); G0' also refuses row 1 | 1.09e-2 m | 0.120 m | 0.875 m | 1.56 m |
| R4 | `gas-order = 1` (cell value sampled) | 1.12e-2 m | 1.37e-2 m | 1.05e-2 m | 1.36e-3 m |

`diam × √2` in the input is refused before G2' ("none of the table diameters"): the oracle reads St
from the diameter. R6, `time-end` absent (St 10): no `snapshot-A.dat`, G0' fails. R7, `out-time =
off`: 10-column rows, G0' refuses ("no time column"). Faces 1–4 retagged 300 is not a RED run: by
QUARTER no parcel of any St comes within h/2 of the boundary (max |x| 1.44, |y| 1.82 m at St 10),
so the mirrored ghosts are never sampled; that mirror is gated by `test_ghost_bc` and
`wall-approach`.
