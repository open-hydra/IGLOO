# swirl-wedge — solid-body swirl on an axisymmetric wedge (W-plan item B)

**Purpose.** The quantitative gate on IGLOO's 2.5D wedge path: a swirling parcel on a
one-cell-thick axisymmetric wedge (`mesh2D` + `axisym`), where the fold (`axisymFold`,
`src/lib/obj_bc.f90`) rotates position AND velocity back into the ±delthe/2 sector at every
segment end. `wedge-fold` (infrastructure) pins the fold *behaviourally* (in-sector, no
stall, exit); this case pins the *numbers* — radius, azimuthal velocity, radial velocity and
the unwrapped azimuth against the exact cylindrical ODE — and measures the 2.5D frame error
the W-plan quantifies (§2) and item E is meant to remove. Registered `swirl-wedge`
(`tests/CMakeLists.txt`, labels `e2e;standard;drag;gas;axisym`).

**Fixture** (`tools/make_wedge_case.py`, all generated; no literature constant). Annular
wedge `x ∈ [0, 2]` (200 cells, dx = 0.01), `r ∈ [0.19, 0.99]` (40 cells, dr = 0.02; cell
centres on round radii), one cell in the azimuth with the k-planes at ∓0.5° (MOSE convention
⇒ `delthe = +0.01745329`). The annulus starts at r = 0.19 so there is no axis, no
`nodeOnAxis` path, and the parcels — which drift outward — never hit the inner face. Gas in
the θ = 0 frame: `U = 1`, `V = 0`, `W = Ω·y_c` with Ω = 0.2 rad/s and `y_c` the vertex-mean
cell centre `r̄_j cos(δ/2)`, so the ord2 bilinear rebuild returns `Ω·y` exactly at any
interior point (`gas_reconstruction` E2); ρ = 1.2, T = 300, μ = 1.8e-5, uniform. `bc.txt`:
f1/f2 400 (outlet), f3/f4 301 (wall), f5/f6 200 (the wedge faces); **no inlet family** —
injection is DB from `input.ini`, and this is the first fixture whose `bc.txt` carries no
401 cell (the reader takes it). The gas ghost fill for 400/301 is `IO.f90`'s default
quadratic extrapolation, exact for this linear field — a non-linear swirl would pick up
extrapolation error at f3/f4. `properties.dat` is vie-plait's with Density 1800 ⇒
`τ = ρ_p d²/(18μ) = 1800·3.6e-7/3.24e-4 = 2.000 s`, St = τΩ = 0.4. Three parcels on one
station (x0 = 0.105, cell 11 centre, θ = 0) with `up = 1 = U` so `u ≡ 1` and
`t = (x − x0)/1` is the clock: P1 co-rotating (r0 = 0.5, w0 = Ωr0 = 0.1), P2 spin-up
(r0 = 0.5, w0 = 0), P3 over-spun (r0 = 0.3, w0 = 2Ωr0 = 0.12). Flight T = 1.895 s;
segment length `Δt_seg = dx/U0 = 0.01 s` (one dual crossing per segment); overshoot per
sector `θ̇Δt_seg/δ ≈ 0.115`.

**Oracle** (`check.py`, derived in-file). In cylindrical coordinates about x, a Stokes
parcel (`a = (g − v)/τ` exactly, Cd = 24/Re) in the true field `g = (U0, −Ωz, +Ωy)` obeys
`ṙ = v_r`, `θ̇ = w/r`, `v̇_r = −v_r/τ + w²/r`, `ẇ = (Ωr − w)/τ − v_r w/r`. Classical RK4 from
the exact injection state, h = T/20000, with two refutations run before every gate: halving
h moves the end state < 1e-12, and RK4 on the 6-state Cartesian system with the true field
agrees with the cylindrical form < 1e-12 (measured 3.7e-15 / 2.9e-15 with DOP853 in the
design study). Measured from `trajectories-A.dat` (`7F12.6,2E13.6E2,I8`, no time column):
`r = hypot(y,z)`, `v_r = (yv + zw)/r`, `w_p = (yw − zv)/r`, `θ_raw = atan2(z,y)`,
`θ_unwrapped = θ_raw + n·delthe` per detected fold (jump < −delthe/2, `n = round`).

**What the code does and the budget the tolerances absorb** (W-plan §2). Between folds
the code samples the θ = 0 Cartesian components `(U0, 0, Ωy)` at the parcel's `(x, y)`
without rotating them to its azimuth: the `−Ωz` component is missing (frame error, first
order in θ) and the fold only re-sectors at segment ends. Mean spurious radial acceleration
`Ω w Δt_seg/(2τ)`, i.e. a relative **outward bias `Δt_seg/(2τ) = 2.5e-3`** on the
centrifugal drift; `w` settles low by `δ²/12 = 2.5e-5`. A segment-wise model of exactly that
behaviour (scratch `swirl_design.py`, DOP853 rtol 1e-12, re-run 2026-09-17) predicts, and
the run reproduces:

| parcel | folds (oracle sectors) | max\|Δr\| model / run | max\|Δw\| | max\|Δθ\| | max\|Δv_r\| |
|---|---|---|---|---|---|
| P1 | 21 (21.09) | 6.8e-5 / **6.75e-5** (end +6.7e-5) | 6.6e-6 / 6.9e-6 | 3.2e-5 / 3.1e-5 | 6.3e-5 / 6.3e-5 |
| P2 | 8 (7.66) | 5.8e-5 / **5.80e-5** (end +5.8e-5) | 1.5e-6 / 1.9e-6 | 1.3e-5 / 1.3e-5 | 6.4e-5 / 6.4e-5 |
| P3 | 33 (32.80) | 7.5e-5 / **7.52e-5** (end +7.5e-5) | 9.5e-6 / 9.8e-6 | 8.1e-5 / 8.2e-5 | 5.8e-5 / 5.8e-5 |

Every residual is positive in r (outward, as predicted) and within 3 % of the model (the
`Δw` differences are the 1e-6 print floor). The gate is therefore a **derived budget**, not
the print floor; once item E (exact 2.5D sampling: rotate the sampled (v, w) to the parcel's
azimuth and sample at (x, r)) lands, these tolerances tighten to the floor and this case is
the one that goes red with E reverted.

**Gates** (`check.py`):

| id | assertion | tolerance |
|---|---|---|
| G0 | 3 parcels; row 1 = `input.ini` state (5e-7); `\|u − 1\| ≤ 1e-6` every row; Tp = 300, d constant; ≥ 100 rows; last x ≥ 1.9 | premise / clock |
| G1 | `\|θ_raw\| ≤ delthe/2 + 1e-6/r` every row (rows are written after the fold) | print floor |
| G2 | folds ≥ 5 per parcel **and** `\|N_folds − round(θ_oracle(T)/delthe)\| ≤ 1` | non-vacuity: Ω = 0 gives 0 folds |
| G3 | `\|r − r_oracle(t)\|` | 2.5e-4 (≈ 3× model) |
| G4 | `\|w_p − w_oracle(t)\|` | 3.0e-5 (≈ 3× model) |
| G5 | `\|θ_unwrapped − θ_oracle(t)\|` | 2.5e-4 (≈ 3× model) |
| G6 | `\|v_r − v_r,oracle(t)\|` | 2.0e-4 (≈ 3× model) |

The diagnostic line per parcel (max residuals, end Δr with sign, fold counts) is printed on
PASS; a first run landing outside ~2× the table above means the fold-cadence model is wrong
and the tolerances must be re-derived before trusting the gate.

**Proven RED (2026-09-17, `b8a1e2f` tree + this fixture).**
- R1 — fold rotating **position only** (`obj_bc.f90` axisymFold, the `stateVar(4:6)` rotation
  commented out, rebuilt): G1/G2 stay green (the parcel is still re-sectored, 22/8/36 folds),
  G3–G6 red at the RED-model magnitudes — P1 `Δr` 2.6e-2 (104× tol), `Δw` 2.2e-3 (75×),
  `Δθ` 1.1e-2 (43×), `Δv_r` 2.4e-2 (120×); P2 2.6e-3 (10×) / 1.3e-4 (4×) / 3.4e-4 (1.4×) /
  5.0e-3 (25×); P3 4.5e-2 (180×) / 6.1e-3 (204×) / 5.1e-2 (205×) / 3.4e-2 (169×). The RED
  signature is inward: a position-only fold injects `−w·delthe` per fold and the parcel
  stops spiralling. Restored and rebuilt: green, trajectories identical to the reference.
- R2 — **no swirl** (`--omega 0`, `wp = 0 0 0`, oracle at Ω = 0): G3–G6 trivially green,
  G2 red on every parcel (`only 0 folds`) ⇒ G2 is the guard that makes the case non-vacuous.

**Floor** (OMP 1 and 5, three runs each): `trajectories-A.dat` identical as a sorted
multiset across all six runs and the gate output byte-identical; only `source.tec` moves at
OMP 5 (the known `!$OMP ATOMIC` re-association floor; not gated). One thread per parcel ⇒
the residuals are deterministic per binary.

**Run signature** (`run_out.txt`): `Block 1 size = 200 40 1`, `Axisymmetric wedge: delthe =
0.01745329 rad`, `2D path: axisymmetric wedge, 2.5D`, `block 1: max|W| = 1.960E-01 -- SWIRL`;
no `stuck in cell` / `no net progress`; 193/192/194 rows, all three exit at x = 2.0.

**Not pinned here.** The Eulerian/source outputs under swirl (W-plan Q7: `euler.tec` carries
no `w_p`, the momentum source lacks the azimuthal component) — `out-file = S` is written and
ignored. The fold cadence for a parcel that does not cross dual cells in x (none here).
