# INFO — evaporation/evap-breakup (ODE model 4 right-hand side at one state, unit)

`test_evap_breakup` calls the production `rhsEvapBreakup` (set up through `setupRHS`, `packAuxVars`, the
drag and heat selectors, no gas interpolation) at one droplet state — d 30 µm, ρ_l 2950, T_p 300 K in a
600 K gas, number rate n = 7.411426900e7 1/s (the `evap-breakup-box` parcel's) — and compares `F(8)` and
`F(9)` with rates typed from the papers, not called from the kernels:

    ṁ_ref = −2π d (k_g/c_pg) ln(1 + B_T),  B_T = c_pg (T_g − T_p)/L_v                  [God53], [Spa53]
    F9_ref = −3 n (d_s − d)/(d τ),  τ = C_s (d/2)/u √(ρ_l/ρ_g),  d_s = σ²/(ρ_g μ_g u³)   [RD87] stripping

Citation tags resolve in [../../REFERENCES.md](../../REFERENCES.md).

## Tests

| id | state | assertion (rel 1e-12) | measured | RED with `F(8) = n·ṁ_evap` |
|---|---|---|---|---|
| U1 | slip 0 (breakup dormant), Y_inf 0 | `F(9) = 0`; `F(8)/ṁ_ref = 1` | 5.8e-16 | ratio **7.41143e7** (= n) |
| U2 | slip 10 m/s (We_r 30, Re 20, stripping), Y_inf 1 (evaporation frozen) | `F(9)/F9_ref = 1`; `F(8)/(−m F9_ref/n) = 1` | 6.6e-16, 5.3e-16 | `F(9)` 6.6e-16 (the breakup rate itself is right); `F(8)` = 0, rel err **1.0** |
| U3 | U2 with Y_inf 0 (both) | stream-mass identity `n F(8) + m F(9) = n ṁ_ref` | 2.7e-15 | rel err **7.41142e7** |
| U4 | U1 with `blowing = LK` against U1 | the convective heat `F(7) m − F(8) L_v/c_p` scales by `f₂ = b/(e^b − 1)`, `b = −(3/2) Pr τ_d ṁ/m`, `τ_d = ρ_l d²/(18 μ_g)` (`[MHB98]` eqs. 17–19, typed in the test), rel 1e-13; `F(8)` unchanged | 8.9e-16 (b 0.919, f₂ 0.610) | — |

U4 is the model-4 call site of the blowing factor, which no end-to-end case reaches (`mhb98-water-blowing` runs
model 2): with that site skipping f₂ U4 fails by 0.64, with b doubled in `blowingFactor` by 0.43.

Every state also asserts `F` finite. With the breakup share dropped from `F(8)` (everything else as is)
U1 passes and U2/U3 fail (1.0 and 18.7): U1 isolates the number-rate factor, U2 the breakup share.

The Reitz-Diwakar arm ignores `acc`, so the drag law is inert here (`Stokes` is selected only because
`interphase` needs a valid selector). A KH-RT state would make `acc` enter through `g_t` and the drag law
would have to be pinned deliberately.
