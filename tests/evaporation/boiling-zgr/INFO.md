# INFO — evaporation/boiling-zgr (boiling of a superheated drop: the Zuo–Gomes–Rutland branch, unit)

`boiling = ZGR` replaces the boiling clamp when the drop's mean temperature is at or above the boiling temperature at
the local gas pressure, `T̄ ≥ T_B(p)`, `p = ρ_g R_g T_g` (Clausius–Clapeyron's inverse, or the inverse of the material's
`Psat` column) `[ZGR00]`, `[OFLEB]`. The surface sits at `T_B` and the vapour leaves there:

    ΔT = T̄ − T_B,  ΔT_f = max(ΔT, 0.5 K)
    α_s = 760 ΔT_f^0.26 (ΔT_f < 5 K), 27 ΔT_f^2.33 (< 25 K), 13800 ΔT_f^0.39   [W m⁻² K⁻¹]
    G_f = α_s ΔT_f π d²/L_v                     (a P2T drop: 10 π k_l d ΔT_f/L_v, the parabolic conduction)
    A = c_p,g (T_g − T_B)/L_v,  T_f = T_B + (T_g − T_B)/3,  k_f = k_g (T_f/T_g)^0.7,  B = π k_f d Nu/c_p,g
    G (1 + G_f/G) = B ln(1 + A (1 + G_f/G))      (G = 0 for A ≤ 0),   S = G + G_f
    Q_s = S c_p,g (T_g − T_B)/(exp(S/B) − 1),   ṁ = −S,   m dh̄/dt = Q_s − S (L_v + h_l(T_B) − h̄)

with `L_v` at `T_B` (the `Lv` column when the material has one) and `Nu` the heat law's. `test_boiling_zgr` calls the
production `zgrBoiling`, `alphaBoil`, `boilingTemperatureCC`, `inverseTableValue`, `tableValue` and, through
`setupRHS`/`packAuxVars`, `rhsEvaporation` (ODE model 2) and `rhsEvapBreakup` (model 4), against values and identities
typed in the test: water (L_v 2.26e6 J/kg, M_v 18, c 4186, ρ 958, d 50 µm) in the box gas (ρ 1.2, R 287, γ 1.4,
k_g 0.026, c_p,g 1004.5), Re 0, Nu 2, CEM below `T_B`. The legs at `A = 0` need `T_B = T_g` exactly: the test steps the
material's `invTboil` by one ulp from its analytic value until `boilingTemperatureCC` gives 600 K (`TsOut` then equals
600 K bitwise).

Citation tags resolve in [../../REFERENCES.md](../../REFERENCES.md).

## Tests

| leg | state | assertion | tolerance | measured | fails on |
|---|---|---|---|---|---|
| (a) | T_g 600 K, T_B 450 K (A = 0.066670, T_f 500 K, k_f 0.0228848), ΔT = 20, 5, 0.5 K | `zgrBoiling`'s −ṁ − G_f against the root of the defining equation typed to ten digits (python, `brentq`): 4.0108994750e-10, 4.6131345055e-10, 4.6190847194e-10 kg/s; the residual `|G (1 + G_f/G) − B ln(1 + A (1 + G_f/G))|/B` | 1e-9 rel.; 1e-13 | 7.3e-12; 2.2e-16 | the blowing factor dropped (G = B ln(1 + A)): 1.517e-1; the film dropped (k_f = k_g): 1.559e-1; 2.33 typed 2.3: 4.20e-1 |
| (b) | T_g 440 K < T_B 450 K (A < 0; T_f 446.67 K) and T_g = T_B (A = 0), ΔT 20 and 0.5 K | ṁ = −G_f (G = 0); `Q_s = G_f c_p,g (T_g − T_B)/(exp(G_f/B) − 1) < 0` at A < 0 (−7.28279057e-5, −8.25402322e-5 W); Q_s = 0 exactly at A = 0 | 1e-14 rel.; exact | 0; 2.5e-16 | Q_s set to 0 below T_B (the reference's): 1.0 |
| (c) | ΔT = 0.2 (under the floor), 1, 10, 30 K | `alphaBoil` against the three laws typed from `[OFLEB]`, ΔT_f = 0.5 K below the floor | 1e-14 rel. | 0 | 2.33 typed 2.3: 6.67e-2 |
| (d) | the adiabatic flash at A = 0 (T_g = T_B = 600 K), T̄ = T_B + 30, + 10, + 1 K, through `rhsEvaporation` | `(F(8)/F(7))/m = c/(L_v − c ΔT)` (dm/dT̄ of the flash with the vapour leaving at T_B: m (L_v − c ΔT) conserved): 1.9611885196e-3, 1.8871667253e-3, 1.8556494463e-3 K⁻¹; `TsOut` = T_B | 1e-13 rel. | 0 | the regime test bypassed (the clamp, CEM): 5.603e-2; the vapour at the mean temperature (no carry): 5.557e-2; the sign of the ln term flipped: 2.127e-1 |
| (e) | p = 0.2, 1, 5 atm; a Psat table of Clausius–Clapeyron on 1-K nodes, 280–500 K | `psat_CC(boilingTemperatureCC(p)) = p` (psat typed); `tableValue(inverseTableValue(p)) = p`, the root inside its 1-K interval | 1e-14; 1e-12 rel. | 6.7e-16; 2.2e-16 | the sign of the ln term flipped: 24 |
| (f) | T̄ = T_B (1 − 1e-12) and T̄ = T_B | below: F(8) equals the clamp run's (boiling = clamp) bitwise; at T_B: F(8) equals `zgrBoiling`'s ṁ bitwise | exact | equal | `>=` typed `>` (the evaporation path at T_B) |
| (g) | (d)'s +30 K state through `rhsEvapBreakup` (Reitz–Diwakar at slip 0: F(9) = 0) | as (d) | 1e-13 rel. | 0 | the regime test dropped from model 4: 5.603e-2; its carry dropped: 5.557e-2 |

`zgrBoiling` solves the defining equation by Brent's method in u = ln(G_f/G) on [−60, 60] to 4ε; a root beyond an end
gives G = B ln(1 + A) (G_f negligible) or G = 0 (G_f dominant).
