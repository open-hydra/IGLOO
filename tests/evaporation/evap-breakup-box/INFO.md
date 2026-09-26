# INFO — evaporation/evap-breakup-box (ODE model 4: breakup, then evaporation to burnout, e2e)

## Reference
- **Breakup:** Reitz, R. D.; Diwakar, R. SAE 870598, 1987 — `[RD87]` (bag D = π, stripping C = 20;
  "the size of the unstable drops was allowed to change continuously … with each change in drop
  size the drop number was changed correspondingly to conserve liquid mass"). The same rule in a
  vaporizing spray: `[Reitz87]` pp. 321–322 (eq. 11, `da/dt = −(a − r)/τ` with `N a³ = N0 a0³`) and
  eq. 21 (vaporization is a drop-radius rate of the drop-number density, so it shrinks drops without
  changing their number).
- **Evaporation:** the stagnant-film d²-law, Godsave `[God53]` / Spalding `[Spa53]`:
  `ṁ = −2π d (k_g/c_pg) ln(1 + B_T)`, `B_T = c_pg (T_g − T_p)/L_v`.
- **Parcel bookkeeping:** PSI-Cell, `[CSS77]` — a trajectory carries a number rate `ṅ_p` and exchanges
  `ṅ_p·Δm` with the cells it crosses.
- No published frozen-gas benchmark combines the two mechanisms (the candidates and why each was
  rejected are in [../LITERATURE_TESTS.md](../LITERATURE_TESTS.md), "Evaporation with ODE breakup"). The
  oracle composes the two kernels, each source-verified elsewhere in the suite (`reitz-diwakar-e2e`,
  `test_breakup_rd`; `d2law`, `test_evap_probes` XE1), along the run's measured `u(x)` and `T_p(x)`.

## The model it verifies
Per droplet of a stream with number rate `n = ṅ_p` (ODE model 4, `rhsEvapBreakup`):

    dm/dt = ṁ_evap − m·(dn/dt)_brk/n        dn/dt = (dn/dt)_brk = −3 n (d_s − d)/(d τ)
    ⇒ d(n m)/dt = n ṁ_evap                  breakup conserves the stream mass
    ⇒ dd/dt = (d_s − d)/τ + 2 ṁ_evap/(ρ_l π d²)

## Fixture
`tools/make_pe_case.py INPUT --we-sweep 30 (×25) --we-convention rad --u-gas 50 --kv 0.8 --sigma 6.0e-5
--lx 0.8 --nx 320 --t-gas 600 --kt 0.5`: a uniform 50 m/s, 600 K gas in a 0.8 × 0.05 × 0.05 m box
(dx = 2.5 mm, 8000 cells); 25 identical 30 µm drops enter at 40 m/s and 300 K (slip 10 m/s, We_r 30,
Re 20: stripping). `phase.txt`/`properties.dat` are `../../common`'s (ρ_l 2950, c_p 1250). σ = 6e-5 N/m
and L_v = 2e5 J/kg are verification inputs chosen so that the breakup window, an evaporation share
inside it and the burnout all fall inside the box; a physical liquid cannot satisfy all three in a
frozen-gas box (drops that break at realistic σ are too large to evaporate within the residence).
`mollify = off`: B2 reads the unsmoothed field.

| quantity | predicted (re-typed closures, before the run) | measured |
|---|---|---|
| d/d0 at the first row past x = 0.04 m | 0.46 | 0.4597 |
| breakup window (rows with We_r ≥ 6.5) | 14–15 rows | ≥ 15 rows |
| burnout (m ≤ 1e-15 kg) | x = 0.543 m, t = 11.15 ms | x = 0.544066 m (all 25) |

## What it verifies

| gate | assertion |
|---|---|
| G0 | witnesses: parents 1..25, injection rows at d0 ± 1 %, u0 ± 1e-3, T0 ± 0.3; Σ outloc column 7 = 25·ṁ_parcel from the case inputs (`krho/(1−krho)·ρ_g U A`) to 1e-6 |
| B1 | closure: Σ `wdot` over `source.tec` = the injected flow to 1e-12 (every drop is consumed in the box) |
| B2 | sign: no cell below −1e-6 of the largest (evaporation is the only mass exchange; breakup moves mass between drops) |
| B3 | end of life: 25 parents end at 0.45 < x < 0.65 m, no `Run_ODESolver err=` line, no `non-finite state` line |
| W1 | breakup acted: the first row past x = 0.04 m has d ≤ 0.55 d0 |
| C1 | both mechanisms: rows with We_r ≥ 6.5, RK4 (40 substeps per row interval) of the composed `dd/dt` along the measured u, T_p (linear between rows, dt = dx/ū), RD87 branch chosen per state, anchored at the injection row; every row within 2e-3; ≥ 20 parcels with ≥ 12 rows |
| C2 | evaporation alone: rows with We_r ≤ 5.5, re-anchored at the first one, the d²-law integral along the measured T_p with the `d2law` case's error budget (E13.6 truncation, Richardson on the T_p sampling, T_p truncation, integrator floor); ≥ 20 parcels with ≥ 5 rows and ≥ 3 % d² loss |

The 5.5–6.5 gap keeps the rate discontinuity at We_r = 6 out of both rate gates.

**Why the mass closure alone cannot be the RED gate:** `computeSource` gives both ends of every segment
from the same state, so Σ `wdot` telescopes to (first in − last out) plus the consumed remnant whatever
the RHS integrated in between; and outloc column 7 of a model-4 parcel is its injected flow, not its
exit flow. The rate gates C1/C2 and the witnesses B3/W1 are the ones that see a wrong mass equation.

## Tolerances
- C1 `2e-3`: ≈ 3× the Richardson term of the T_p/u sampling (6.3e-4 at the window end, derived on the
  exact solution sampled at every cell crossing and printed F12.6/E13.6). Measured worst **1.39e-4**.
- C2: the d2law case's budget, unchanged. Measured worst resid/tol **0.512**.
- B1 `1e-12` against the case-input flow (no printed value enters): measured **5.4e-16**.
- B2 measured: min `wdot` = 0.0 (no cell below zero), max 6.16e-5 kg/s.

## Falsification (measured)
| solver | outcome |
|---|---|
| unfixed (`F(8) = n·ṁ_evap`, no breakup share) | every parcel lost at injection: SDIRK4 gives up at t = 1.88e-10 s (predicted 1.9e-10 s), 25 `Run_ODESolver err=-1` lines; B1 Σ `wdot` = 0, resid **1.000**; B2 no deposit; B3 0/25 in band; W1 0/25; C1, C2 0 parcels |
| breakup share dropped (`F(8) = ṁ_evap`) | B1 Σ `wdot` = −1.112 kg/s (the gas loses mass to the breakup cells), B2 min/max **−8.5**, W1 d/d0 = 0.966, B3 parcels leave alive at x = 0.8, C1 **0.881** (predicted 0.68 on a shorter window) |
| model 4 left out of the burnout test | B3 only: 25 `Run_ODESolver err=-1` lines at t = 12.3 ms, x = 0.6025 m; B1 still closes (the failure exit hands the remnant to the gas) |
| fixed | all seven PASS |

## Pass criterion
All seven gates PASS. The witness bands (G0 anchors, B3's x band, W1's 0.55) are centred on the measured
run, not on a model: a burnout moving out of the band is a finding about the solver, not a reason to widen it.
