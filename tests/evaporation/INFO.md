# INFO — evaporation (family C)

Verifies the evaporation mass ODE (and coupled heat) under the production
closures: d²-law=1, CEM=2, CEM-B=3, ASM=4, TC=5. All property-freeze
flags are forced `.false.` ⇒ constant
gas/liquid properties.

Citation tags resolve in [../REFERENCES.md](../REFERENCES.md).

## Tests

| id | mode | source | doi | locus | compared | settings_provenance | status |
|---|---|---|---|---|---|---|---|
| C0 | analytic (manufactured) | derived | n/a | "constant-area-flux" target — **not directly realizable** through d2law (gives `mdot ∝ d`, not `d²`). Implemented as a manufactured/forcing note rather than an exercised test | documentation row only; no executable test | n/a | no test (see note) |
| C1 | analytic | `[Lef89]` (d²-law derivation) | ISBN 0-89116-697-5 | LITERATURE `mdot = −2π d k_g/c_pg · ln(1+B_T)` (Nu=2 stagnant film) with frozen `T_p` ⇒ `d(d²)/dt = −K`, `K = 8 k_g ln(1+B_T)/(ρ_l c_pg)`. **NB:** the oracle uses the literature 8, so it fails a half-rate code (`−π d …`, `K=4`) | linear regression of `d²(t)` vs t: slope match to K, `R² > 0.999`, lifetime `t_life = d₀²/K` match | frozen `T_p`, constant `c_pg, k_g, ρ_l, L_v`; `B_T = c_pg(T_g−T_p)/L_v` from inputs | e2e: `d2law`, `d2law-line` |
| C2 | analytic | `[Lef89]` (wet-bulb / heat-up theory) | ISBN 0-89116-697-5 | heat-up phase → asymptotic wet-bulb `T_wb`; `T_wb` independent of `T_p₀` and `d₀` for a given (gas state, fuel) | terminal `T_p` matches independent `T_wb` from psychrometric/iteration; invariance across (T_p₀, d₀) sweep | constant gas state, single fuel, no breakup | e2e: `mhb98-water` (wet-bulb identity), `lk-neq` (heat-up) |
| C3 | analytic + hybrid | derived (Abramzon–Sirignano / ASM family) | see note | ASM `F(B) = (1+B)^0.7 · ln(1+B)/B`; Sh*/Nu* modified; ψ–B_T iteration. Limits: `B_T → 0` collapses to d²-law (C1) | C1-limit collapse within tol; `relL2 < tol` vs an independent oracle (self-contained RK4 of ASM mass-energy system); `p_obs ≈ 4 ± 0.3` for L2 | constant properties, fuel/gas state ranges spanning low/moderate B_T | e2e: `asm-box` (d² and T along the measured state); unit: EV10, EV11 (`test_evaporation`), XE3 (`test_evap_probes`); the C1-limit collapse and the order check: no test |
| C4 | analytic | derived (energy & mass balance) | n/a | global energy balance `Q_in = m_v L_v + m c_p (T_wb − T_p₀)` (or equivalent); global mass balance `m(t) + m_v(t) = m₀` | both balances close to ≤ 10·ε_mach (relative) over integrated trajectory | constant properties; integrator runs to `d → ~0` | partial: mass telescoping in `d2law/check.py`; the energy half: no test |

## Notes

- **C0:** a known modelling mismatch (d²-law's
  `mdot ∝ d` precludes a clean constant-area-flux mode): no executable test.
- **C3 source:** the ASM/A–S family combines several papers (Abramzon &
  Sirignano 1989 for the core formulation, `[AS89]`; Sazhin reviews for compact summaries).
  `asm-box` types the rate law from `[MHB98]`'s statement of the model and the gas-side heat
  from `[AS89]` eq. 20.
- **C1 frozen-`T_p`:** the d²-law closed form assumes `T_p` constant;
  `d2law-line` realizes it with a particle `cp` 100× and a density 2.7× the `d2law` table's
  (`ρ_p·cp_p` 271×; `T_p` drifts 0.77 K over the residence), and `d2law` integrates the law
  along the measured `T_p`.

## Unit tests

- `test_evaporation.f90` (ctest `test_evaporation`, GREEN): function-level
  checks of `evaporation()` called DIRECTLY with correctly-ordered gas args.
  EV1 CEM Re=0 vs independent Spalding chain (psat_CC→Xs→Ys→BM, Sh=2), 2e-16;
  EV2 Ranz-Marshall Sh ratio exact; EV3 CEM-B 1/3-rule film chain at Re=0,
  6e-19. Confirms the models at function level (the production call site is
  XE2's).
  EV10 ASM near equilibrium (surface state Y_s = 0.25 exactly, B_M = k·2⁻⁵⁵/0.75 from
  3.7e-17 to 1.04e-2, three φ/Re settings): ṁ and Q_G against the same algorithm in
  quadruple precision from the same inputs, 1e-12; EV11 the corner B_M ≤ 1.1e-16 (and
  φ = 0.3 at k = 4, 8): Q_G against the Stefan-free limit π·d·k_g·Nu₀·(T_g − T_p), 1e-12.
- `test_evap_probes.f90` (ctest `test_evap_probes`, ordinary gate — transcription pin; a failing
  probe turns it red): XE1 d2law/literature ratio = 1 within 1e-6 (a d2law that codes half the
  rate gives exactly 0.5); XE2 single REAL `rhsEvaporation` evaluation (production call site, gas
  packed per contract) → F(8)=mdot < 0 (a scrambled call site gives F(8)=mdot≡0); XE3 ASM
  gas-side heat Q_G > 0 with the heat override set and mdot < 0 (a Q_G = mdot·(Lv + cpv·dT/BT),
  the latent heat counted twice, is negative).
- `test_vapour_tables.f90` (ctest `test_vapour_tables`): the `Lv` and `Cpv` columns of
  `properties.dat` at the droplet temperature. One production right-hand side per leg at
  T_p = 337.5 K with one property tabulated on 280–400 K and its packed scalar poisoned (1e30),
  against the same evaluation with the table off and the scalar set to the table's value at
  T_p, equal to 1e-15 relative (measured 0) and the drop evaporating: VT1 CEM, the latent sink
  F(7); VT2 d²-law, the rate F(8) (B_T); VT3 TC, `Cpv` in F(7) and F(8); VT4 CEM + Reitz-Diwakar
  (ODE model 4), the latent sink of ṁ_evap. RED on one-line mutations of `Lib_RHS.f90`: the
  model-2 `Lv` override deleted (VT1 1.9e23, VT2 F(8) = 0), the model-4 one (VT4), the `Cpv`
  override (VT3 2.4 and 0.67), the model-2 sink back to the packed scalar (VT1 only, VT2 green)
  and the model-4 one (VT4 only).

## Evaporation with ODE breakup (ODE model 4)

- `evap-breakup/test_evap_breakup.f90` (unit): one production `rhsEvapBreakup` call per state against the
  d²-law and Reitz-Diwakar rates typed from the papers — see [evap-breakup/INFO.md](evap-breakup/INFO.md).
- `evap-breakup-box` (e2e): breakup then evaporation to burnout, the composed-rate oracle and the source
  closure — see [evap-breakup-box/INFO.md](evap-breakup-box/INFO.md).
- `d2law-brk-dormant` (e2e): the d2law case under model 4 with the breakup rate exactly zero.

## Heat conduction inside the drop (`liquid-conduction = P2T`)

- `liquid-p2t/test_liquid_p2t.f90` (unit): the parabolic temperature profile through the right-hand sides of
  models 2 and 4 — the closed form at Re = 0, the surface balance with the carry of the evaporated mass, χ(Pe_L),
  the breakup viscosity as the liquid one, radiation inside the balance — see [liquid-p2t/INFO.md](liquid-p2t/INFO.md).
- `temp-relax-p2t` (e2e): temp-relax's drop with no mass transfer, τ = 1.2 τ_ITC and the trajectory's
  surface-temperature column — see [temp-relax-p2t/INFO.md](temp-relax-p2t/INFO.md).

## Boiling of a superheated drop (`boiling = ZGR`)

- `boiling-zgr/test_boiling_zgr.f90` (unit): the boiling branch, its heat-transfer coefficient, the boiling
  temperature at the pressure and the adiabatic flash through models 2 and 4 — see [boiling-zgr/INFO.md](boiling-zgr/INFO.md).
- `boiling-zgr-flash`, `boiling-zgr-hot`, `boiling-zgr-switch` (e2e): the flash at A = 0, the boiling in a hotter gas,
  the switch into evaporation at the boiling point and the burnout — see their `INFO.md`.
