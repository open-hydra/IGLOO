# INFO — breakup_tab

Verifies the TAB (Taylor Analogy Breakup) model. **Not SDIRK-integrated** —
`y(t)` is advanced by closed-form damped-oscillator update `YupdateTAB`
in `solout` (event-based). This family is direct-call verification of
`breakupEvent`, not integrator-driven.

Citation tags resolve in [../../REFERENCES.md](../../REFERENCES.md).

## Tests

| id | mode | source | doi | locus | compared | settings_provenance | status |
|---|---|---|---|---|---|---|---|
| D1 | analytic | `[ORA87]` O'Rourke & Amsden 1987 | DOI:10.4271/872089 | TAB coefficients: `rdt = ½·Cmu·μ_l/(ρ_l r²)`, `ω² = Comega·σ/(ρ_l r³) − rdt²`, `WeCr = (r ρ_g v²/σ)/WeCrit`; default `bp = [Comega=8, Cmu=5, WeCrit=12, nSpread]`; O'Rourke–Amsden mapping `C_k=Comega`, `C_d=Cmu`, `WeCrit ≈ C_b C_k / C_F` | (a) `y(t)`, `ẏ(t)` vs independent reference damped-oscillator integrator (in `verif_oracle.f90`, **not** calling production `YupdateTAB`): `relL∞ < tol`; (b) trigger time of `y > 1` and asymptotic `y_∞` correct; (c) the product radius drawn from the size law's mass distribution: its moments and the Sauter radius the parcels represent, from a fixed-seed stream; (d) subcritical We (`We < WeCrit`) ⇒ no breakup event fires | constant gas state spanning sub- / mid- / supercritical We; the production `bp` at the [ORA87] defaults set in the test driver (the INI→`bp` reader: [ini_pipeline](../../infrastructure/ini_pipeline/INFO.md) IP1) | **GREEN** (`test_breakup_tab.f90`, ctest `test_breakup_tab`): (a) done twice — code vs analytic 3.6e-15 AND RK4-of-raw-ODE vs analytic 3.8e-15 (triangulation); (b) tb via event-flag bisection vs RK4 first `y=1` crossing, err 7.3e-10 s inside the damping-neglect bound; (d) rest-drop sweep flips exactly at `We_r=6=bp(3)/2`. (c) product-size moments: in `test_tab_moments.f90` (ctest `test_tab_moments`) — 20000 seeded from-rest events (TABM), sample E[r]/E[r²] vs the moments m4/m3, m5/m3 of an independent truncated-Rosin-Rammler number pdf, i.e. its mass distribution (own Simpson quadrature + bisection on the D32=rs condition; k=2.8333 energy balance re-derived), 5·CLT-se + fixed-point-stop tolerance: E[r] off by 2.0e-8 against 3.0e-7, E[r²] by 6.1e-13 against 2.3e-11; TABM-S the represented Sauter radius N/Σ(1/r) = 1.0018 r_s against the window's D32 = r_s, 1.8e-3 against 9.8e-3 relative; also gates child∈[2rMin,dp) and npdot·(d_old/d_new)³ per event. RED with the radius drawn from the number pdf: E[r] off by 7.6e-6 (25× the tolerance), the represented Sauter radius 0.737 r_s (27×). INI→bp pipeline: in [../../infrastructure/ini_pipeline](../../infrastructure/ini_pipeline/INFO.md) (IP1 pins bp=[8,5,12,3.5], bpMethod=2, bpScale through the production reader) |
| TABM1 | analytic + stochastic | `[KIVA2]` (eqs. (30)–(31), the exponential product law, r32 = 3 r̄; App. D and App. F, the radius drawn from the mass distribution) | LA-11560-MS | TAB resize method 1 (`Lib_Breakup.f90` TABmodel): the size law is chi-square with ν = 2 degrees of freedom times rSize = r_s/(ν + 4) = r_s/6 (KIVA-II's exponential with r̄ = r_s/3; untruncated Sauter radius r_s: E[X³]/E[X²] = ν + 4), the radius drawn from its mass distribution, chi-square with ν + 6 = 8, rejected into [0.01 r_s, r) | sample E[r], E[r²] of 20 000 from-rest events through `breakupEvent` (`test_tab_moments.f90`, after the method-2 leg, same stream) against the mass-pdf moments m4/m3, m5/m3 of the χ²_2 pdf x^(ν/2−1) e^(−x/2)/(2^(ν/2) Γ(ν/2)) scaled by rSize and truncated to the same window, by the test's Simpson quadrature; 5 CLT standard errors + 1.5e-4 / 3e-4 relative; TABM1-S the represented Sauter radius N/Σ(1/r) against the window's number-pdf D32, 5 CLT standard errors + 1.5e-4 relative; every child in the window, npdot·(d_old/d_new)³ per event | k = 2.8333, rSize = 5.882e-6 m, mass-pdf E[r] = 4.491e-5 m, window D32 = 0.97898 r_s | **GREEN**: E[r] off by 3.1e-7 against 7.2e-7, E[r²] by 2.6e-11 against 7.4e-11; TABM1-S 0.98996 r_s against 0.97898 r_s, 1.1e-2 against 2.5e-2 relative, 2.3 standard errors from the oracle on this fixed stream (the leg closest to its tolerance). RED with the number draw at ν = nint(k) = 3: E[r] off by 3.0e-5 (41× the tolerance), the represented Sauter radius 0.180 r_s (33×) |

## Notes

- Coefficient mapping `[ORA87] ↔ IGLOO`:
  - `[ORA87] C_k = 8` ↔ IGLOO `Comega = 8`
  - `[ORA87] C_d = 5` ↔ IGLOO `Cmu = 5`
  - `[ORA87] WeCrit ≈ 12` ↔ IGLOO `WeCrit = 12`
- The formulation is verified here, the value pipeline from the INI to `bp` by the
  production reader in [ini_pipeline](../../infrastructure/ini_pipeline/INFO.md) (IP1).
- Product sizes are drawn from the mass distribution of the size law (chi-square for
  method 1, Rosin–Rammler through `GammaS`, `[MT00]`, for method 2): a parcel keeps its
  mass flow, so it represents the law. `test_tab_moments` draws from its own fixed-seed
  stream (`tRng`, passed to `breakupEvent`) for reproducibility. Moments (not
  realisations) are compared to closed-form values.
