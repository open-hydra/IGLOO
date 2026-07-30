# INFO — breakup_tab (family D, test D1)

Verifies the TAB (Taylor Analogy Breakup) model. **Not SDIRK-integrated** —
`y(t)` is advanced by closed-form damped-oscillator update `YupdateTAB`
in `solout` (event-based). This family is direct-call verification of
`breakupEvent`, not integrator-driven.

Schema: see verification-phase1-plan.md §2.5.
Citation tags resolve in [../../REFERENCES.md](../../REFERENCES.md).

## Tests

| id | mode | source | doi | locus | compared | settings_provenance | status |
|---|---|---|---|---|---|---|---|
| D1 | analytic | `[ORA87]` O'Rourke & Amsden 1987 | DOI:10.4271/872089 | TAB coefficients: `rdt = ½·Cmu·μ_l/(ρ_l r²)`, `ω² = Comega·σ/(ρ_l r³) − rdt²`, `WeCr = (r ρ_g v²/σ)/WeCrit`; default `bp = [Comega=8, Cmu=5, WeCrit=12, nSpread]`; O'Rourke–Amsden mapping `C_k=Comega`, `C_d=Cmu`, `WeCrit ≈ C_b C_k / C_F` | (a) `y(t)`, `ẏ(t)` vs independent reference damped-oscillator integrator (in `verif_oracle.f90`, **not** calling production `YupdateTAB`): `relL∞ < tol`; (b) trigger time of `y > 1` and asymptotic `y_∞` correct; (c) child-size distribution moments (mean, std) match hand-computed values with `rng_seed=42`; (d) subcritical We (`We < WeCrit`) ⇒ no breakup event fires | constant gas state spanning sub- / mid- / supercritical We; INI `bp` read in test driver and compared against [ORA87] defaults | **partial — GREEN** (`test_breakup_tab.f90`, ctest `breakup_tab`): (a) done twice — code vs analytic 3.6e-15 AND RK4-of-raw-ODE vs analytic 3.8e-15 (triangulation); (b) tb via event-flag bisection vs RK4 first `y=1` crossing, err 7.3e-10 s inside the damping-neglect bound; (d) rest-drop sweep flips exactly at `We_r=6=bp(3)/2`. (c) child-size moments: **DONE 2026-07-02** in `test_tab_moments.f90` (ctest `test_tab_moments`) — 20000 seeded from-rest events, sample E[r]/E[r²] vs an independent truncated-Rosin-Rammler oracle (own Simpson quadrature + bisection on the D32=rs condition; k=2.8333 energy balance re-derived), 5·CLT-se + fixed-point-stop tolerance, both within ~2σ; also gates child∈[2rMin,dp) and npdot·(d_old/d_new)³ per event. INI→bp pipeline: **DONE** in [../../infrastructure/ini_pipeline](../../infrastructure/ini_pipeline/INFO.md) (IP1 pins bp=[8,5,12,3.5], bpMethod=2, bpScale through the production reader) |

## Notes

- Coefficient mapping `[ORA87] ↔ IGLOO`:
  - `[ORA87] C_k = 8` ↔ IGLOO `Comega = 8`
  - `[ORA87] C_d = 5` ↔ IGLOO `Cmu = 5`
  - `[ORA87] WeCrit ≈ 12` ↔ IGLOO `WeCrit = 12`
- The test **reads** `bp` from the production INI parser to confirm the
  defaults — i.e. it verifies both the formulation **and** the value-pipeline
  from INI to `YupdateTAB`.
- Child size sampling is stochastic (ChiSquare / Rosin–Rammler); test fixes
  `rng_seed = 42` for reproducibility. Moments (not realisations) are compared
  to closed-form values.
