# INFO — evaporation/burnout-ord2-slab (an evaporating drop burns out at gas order 2 on a planar slab, e2e)

## What it guards
At gas order 2 the right-hand side samples the gas at the trial position the ODE solver hands it: on a planar
one-layer mesh, `sampleGas2D` → `interp2ndOrder2D`, the bilinear map of the parcel's dual cell. SDIRK4's Newton
iterations may try states far from the accepted one: here, with `[IGLOO-ODE] mass-tol = off`, the mass slot keeps
the default `absolute-tol` 1e-10, above a 4.2e-12 kg drop's mass, the last part of its life is outside the solver's
error control, and the trial states reach positions of order 1e16 m with a negative mass. An unclamped bilinear map
evaluated there extrapolates four equal node values into γ and R of order 1e20 with either sign, and the
evaporation model's guard on the gas heat capacity (`cpg <= 0`) stops the run on a state the solver would have
rejected. Two defences keep the run going: the bilinear sample clamps its local coordinates to the dual cell, and a
right-hand side handed a trial state with no mass left (or a non-finite slot) at a gas of `cpg <= 0` returns the
1e30 penalty that makes SDIRK4 reject the step. The case gates that no trial state stops the run and that every
drop ends on a burnout record with its whole mass handed to the gas.

The default, `mass-tol = burnout`, gives the mass slot `min(absolute-tol, 1e-3 m_b)` = 1e-23 kg here: the trial
states stay in the domain. A sibling run of the same case without the key (`ref-burnout/`) witnesses that path.

At gas order 1 the gas is the cell's own value, the same trial states reach the right-hand side and the run ends
clean: the stop is the gas-order-2 sample's.

## Case construction
`make_fixture.py` writes `INPUT/solfile.tec` and `INPUT/bc.txt` with `tools/make_box_case.py`: a one-layer box,
20 × 20 × 1 cells over 0.1 × 0.1 × 0.005 m (IGLOO's planar 2D path), every face 400, a uniform nitrogen-like gas at
rest: T 350 K, γ 1.399, R 296.8 J/(kg K), ρ 0.9754 kg/m³ (101325 Pa), μ 2.0e-5 Pa s, k 0.0293 W/(m K).
`INPUT/properties.dat` is `mhb98-water`'s (ρ_p 997, c_p 4184, constant; symlinked), `INPUT/phase.txt` the common one.

Four drops by `[IGLOO-BC]` coordinates at x = 0.02 m, y = 0.03, 0.045, 0.055, 0.07 m, z = 0.0025 m (four dual
cells, none on a node line): d₀ = 20 µm, T₀ = 320 K, u₀ = 0.5 m/s along x, 1e-6 kg/s each. Stokes drag stops them
within u₀ τ₀ = 0.554 mm (τ₀ = ρ_p d₀²/(18 μ_g) = 1.108 ms); CEM evaporation with the Langmuir-Knudsen interface, the
water constants of `mhb98-water` (L_v 2.462478e6 J/kg, M_v 18.015, T_boil 373.15 K, Clausius-Clapeyron p_sat), then
burns them out at m_b = min(1e-6 m₀, 1e-20 kg) = 1e-20 kg (m₀ = 4.176e-12 kg, d(m_b) = 26.76 nm). `out-time = on`
puts `t` on every exit record; `out-file = S` writes the source field. `[IGLOO-ODE] mass-tol = off`; `check.py` writes
the sibling deck, the same lines without the `[IGLOO-ODE]` section, into `ref-burnout/` and runs it there.

## What it verifies (`check.py`, inputs, logs and the printed records only)

| id | gate | tolerance |
|---|---|---|
| P0 | no line of `run_out.txt` or `run_err.txt` carries an error stop, an `ERROR`, a non-finite state, a solver failure or the inner-loop cap (the endings that are not a burnout) | exact |
| P1 | every stop leaves an exit record: one per parcel ID 1–4, `bc` 0 (consumed), `t > 0` | exact |
| P2 | the drop stops where drag stops it: `0 < X − x₀ ≤ u₀ τ₀ + 5e-7 m`, Y and Z at the injection point (the drop shrinks, so its relaxation time only falls); a trial state accepted outside the box fails it | half the F12.6 print |
| P3 | the burnout state: `d_p ≤ d(m_b)(1 + 1e-5)`, `d(m) = (6m/(π ρ_p))^{1/3}`; a drop removed before burnout fails it | twice the E13.6 print of `d_p` |
| P4 | the mass reaches the gas: Σ `wdot` over `source.tec` against Σ `mdot` of the exit records | 1e-5 (`tc-hexadecane`'s closure) |
| P0–P4 (ref-burnout) | the same five gates on the sibling run (default `mass-tol`) | as above |
| P5 | the mass slot's tolerance the solver prints: 1e-10 with `mass-tol = off` (this run); `min(1e-10, 1e-3 m_b)` = 1e-23 kg and its burnout mass `m_b` = 1e-20 kg in `ref-burnout/` (exit 0) | 5e-4 (the ES10.3 print) |

Measured: P2 `X − x₀` 5.40e-4 m on every drop; P3 `d_p/d(m_b)` 0.86055, 0.90273, 0.95179, 0.92817 (this run) and
0.87803, 0.79840, 0.87774, 0.87733 (`ref-burnout`); P4 residual 0 in both; the outputs are the same bytes at 1 and
5 threads.

## RED (`check.py` on each)

| solver | what fails |
|---|---|
| the bilinear sample unclamped, no reject of a trial state | the run stops in the first parcel's tail: `evaporation: cpg<=0 (gas args scrambled?)` on stderr, exit 128, no exit record (P0, P1, P2, P3 and P4 fail on the partial output) |
| the same with the mass slot at `min(atol, 1e-3 m_b)` | the same stop with `mass-tol = off`; `ref-burnout/` passes P0–P4 |
| the mass slot at `absolute-tol` whatever `mass-tol` says (no `mass-slot` line) | P5 |

Either defence alone passes P0–P4: the clamp alone (P3 0.86055–0.95179, the numbers above) and the penalty alone
(P3 0.63133, 0.81721, 0.80986, 0.88800). The same case at gas order 1 or with `absolute-tol = 1e-25` passes P0–P4 on
the solver with neither: P3 `d_p/d(m_b)` 0.708–0.838 (gas order 1) and 0.913–0.975 (`absolute-tol` 1e-25).
