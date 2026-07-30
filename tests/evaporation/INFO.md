# INFO — evaporation (family C)

Verifies the evaporation mass ODE (and coupled heat) under the four production
closures: d²-law=1, CEM=2, CEM-B=3, ASM=4, LEB=5. Per plan OQ1 (resolved
2026-06-13), all property-freeze flags are forced `.false.` ⇒ constant
gas/liquid properties.

Schema: see verification-phase1-plan.md §2.5.
Citation tags resolve in [../REFERENCES.md](../REFERENCES.md).

## Tests

| id | mode | source | doi | locus | compared | settings_provenance | status |
|---|---|---|---|---|---|---|---|
| C0 | analytic (manufactured) | derived | n/a | "constant-area-flux" target — **not directly realizable** through d2law (gives `mdot ∝ d`, not `d²`). Implemented as a manufactured/forcing note rather than an exercised test | documentation row only; no executable test | n/a | xfail (deferred, see note) |
| C1 | analytic | `[Lef89]` (d²-law derivation) | ISBN 0-89116-697-5 | LITERATURE `mdot = −2π d k_g/c_pg · ln(1+B_T)` (Nu=2 stagnant film) with frozen `T_p` ⇒ `d(d²)/dt = −K`, `K = 8 k_g ln(1+B_T)/(ρ_l c_pg)`. **NB (FINDINGS D-EVAP-1):** production `d2law` codes half this (`−π d …`, `K=4`); the oracle must use the literature 8 so it fails the half-rate code | linear regression of `d²(t)` vs t: slope match to K, `R² > 0.999`, lifetime `t_life = d₀²/K` match | frozen `T_p`, constant `c_pg, k_g, ρ_l, L_v`; `B_T = c_pg(T_g−T_p)/L_v` from inputs | not started |
| C2 | analytic | `[Lef89]` (wet-bulb / heat-up theory) | ISBN 0-89116-697-5 | heat-up phase → asymptotic wet-bulb `T_wb`; `T_wb` independent of `T_p₀` and `d₀` for a given (gas state, fuel) | terminal `T_p` matches independent `T_wb` from psychrometric/iteration; invariance across (T_p₀, d₀) sweep | constant gas state, single fuel, no breakup | not started |
| C3 | analytic + hybrid | derived (Abramzon–Sirignano / ASM family) | TBD (see note) | ASM `F(B) = (1+B)^0.7 · ln(1+B)/B`; Sh*/Nu* modified; ψ–B_T iteration. Limits: `B_T → 0` collapses to d²-law (C1) | C1-limit collapse within tol; `relL2 < tol` vs Tier-2 oracle (self-contained RK4 of ASM mass-energy system); `p_obs ≈ 4 ± 0.3` for L2 | constant properties, fuel/gas state ranges spanning low/moderate B_T | not started |
| C4 | analytic | derived (energy & mass balance) | n/a | global energy balance `Q_in = m_v L_v + m c_p (T_wb − T_p₀)` (or equivalent); global mass balance `m(t) + m_v(t) = m₀` | both balances close to ≤ 10·ε_mach (relative) over integrated trajectory | constant properties; integrator runs to `d → ~0` | not started |

## Notes

- **C0 deferred:** plan §3 surfaces this as a known modelling mismatch (d²-law's
  `mdot ∝ d` precludes a clean constant-area-flux mode). Logged as `xfail` with
  worth recording on first run.
- **C3 source TBD:** the ASM/A–S family combines several papers (Abramzon &
  Sirignano 1989 for the core formulation, Sazhin reviews for compact summaries).
  Implementer to pick the primary citation at C3 implementation time and add to
  `REFERENCES.md` (suggested tag: `[AS89]`, DOI:10.1016/0017-9310(89)90064-7).
- **C1 frozen-`T_p`:** the d²-law closed form assumes `T_p` constant; per OQ1
  this is realized by holding `T_p` at the wet-bulb value (or any constant
  consistent with the fuel state).

## Unit tests (2026-07-02)

- `test_evaporation.f90` (ctest `test_evaporation`, GREEN): function-level
  checks of `evaporation()` called DIRECTLY with correctly-ordered gas args.
  EV1 CEM Re=0 vs independent Spalding chain (psat_CC→Xs→Ys→BM, Sh=2), 2e-16;
  EV2 Ranz-Marshall Sh ratio exact; EV3 CEM-B 1/3-rule film chain at Re=0,
  6e-19. Confirms the models are sound at function level — the breakage is at
  the call sites.
- `test_evap_probes.f90` (ctest `test_evap_probes`, WILL_FAIL xfail): XE1
  d2law/literature ratio measured = exactly 0.5 (bug A3); XE2 single REAL
  `rhsEvaporation` evaluation (production call site, gas packed per contract)
  → F(8)=mdot≡0 (A4). Suite flips RED when both are fixed → promote + unblock
  C1-C4 and tests/evaporation/d2law.
