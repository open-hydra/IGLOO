# INFO — evaporation/liquid-p2t (heat conduction inside the drop: the parabolic temperature profile, unit)

`liquid-conduction = P2T` takes the liquid temperature inside the drop parabolic, `T(r) = T_c + (T_s − T_c)(r/R)²`
`[DS03]`. The ODE state stays the mean temperature `T̄ = (2T_c + 3T_s)/5`; the heat the liquid absorbs at the surface,
`Q = π d² · 2k_eff (T_s − T_c)/R`, fixes the surface temperature in every right-hand-side call:

    T_s − T̄ = Q_net(T_s)/(10 π k_eff d),   Q_net = Q_g(T_s) + Q_r(T_s) + ṁ(T_s) L_v(T_s)
    k_eff = χ k_l,  χ = 1.86 + 0.86 tanh[2.245 log10(Pe_L/30)]                 [AS89] eq. 39
    Pe_L = U_s d ρ_l c_l/k_l,  U_s = (12.69/32) |ΔU| (μ_g/μ_l) Re^(1/3)/(1 + B_M)  [AS89] eqs. 25, 30, 32

with `χ` once per call at `T̄` (`B_M` of the gas-side model at `T̄`; `χ = 1` without circulation, `Re = 0`). The
gas-side models, the convective heat (and the drag's temperature ratio), the radiative heat and the vapour data take
`T_s`; the evaporated mass leaves at `T_s`, so the mean temperature carries the receding surface's enthalpy:

    m dh̄/dt = Q_g(T_s) + Q_r(T_s) + ṁ (L_v(T_s) + h_l(T_s) − h̄)        (at constant c: m c dT̄/dt = … + ṁ c (T_s − T̄))

`test_liquid_p2t` calls the production `setupRHS`, `packAuxVars`, `rhsEvaporation` (ODE model 2), `rhsEvapBreakup`
(model 4), `p2tSurface` and `chiETC` (no executable, no gas interpolation) and compares with closed forms typed in the
test. The gas is 1.2 kg/m³, μ_g 1.8e-5 Pa s, k_g 0.026 W/m/K, γ 1.4, R 287; the vapour water-like (L_v 2.26e6 J/kg,
M_v 18, T_boil 373.15 K, c_p,v 1900, Le 1); drag Stokes, heat Ranz–Marshall.

Citation tags resolve in [../../REFERENCES.md](../../REFERENCES.md).

## Tests

| leg | state | assertion | tolerance | measured | fails on |
|---|---|---|---|---|---|
| (a) | Re 0 (slip 0), Y_inf 1 (ṁ 0), T̄ 300 K, T_g 600 K, d 11.89 µm, ρ 2950, c 1250, k_l = k_g (ζ = Nu k_g/(2k_l) = 1) | `F(7) = π d k_g Nu · 5(T_g − T̄)/((5 + ζ) m c)` | 1e-13 rel. | 2.2e-16 | a right-hand side that ignores `liquid-conduction`: F(7) is 6/5 of the closed form (rel. error 0.2000); the surface temperature replaced by T̄ after the solve: 0.2000 |
| (a′) | as (a) | `p2tSurface` and `TsOut` give `T_s = (5T̄ + ζT_g)/(5 + ζ) = 350 K` | 1e-13 rel. | 8.7e-19 | a solve that returns T̄: 1.429e-1 |
| (b) | (a) with k_l = 1e12 | `F(1:8)` of P2T equal those of ITC (`T_s − T̄` = 1.5e-12 K) | 1e-12 rel. | 4.7e-15 | — (a control) |
| (c) | `chiETC` at Pe_L = 0, 30, 300, 1e12 | 1, 1.86, 1.86 + 0.86 tanh(2.245), 2.72 (eq. 39 typed in the test) | 1e-14 abs. | 4.4e-16 | the constant 2.245 typed 2.225: 7.70e-4 |
| (d) | CEM, T̄ 320 K, T_g 800 K, slip 6 m/s (Re 20), d 50 µm, ρ 1000, c 4180, k_l 0.6, μ_l 1e-3; the temperature state, then the enthalpy state (h = cT) | control: `|G(T_s)| ≤ 4ε T̄` with `K` formed in the test, `T̄ < T_s < T_g`; gate: `F(7) = (Q_g + ṁ (L_v + c (T_s − T̄)))/(m c)` (h state: `/m`), `Q_g` and `ṁ` from `interphase` and `evaporation` at `TsOut` | 1e-13 rel. | `T_s` 323.682316697 K, `G` 1.1e-14 K (0.04 of 4ε T̄); F(7) 2.640358131040e4 K/s (2.2e-16); h state 1.103669698775e8 J/kg/s (2.2e-16) | the carry dropped: 1.210e-3; the carry with T̄ for T_s: 1.210e-3 |
| (e) | T_g ∈ {300, 600, 1500} K × T̄ ∈ {280, 350, T_boil − 0.1, T_boil + 5} K × d ∈ {1, 100} µm × k_l ∈ {0.01, 0.6}, slip 6 m/s; d²-law, CEM, CEM-B, ASM, TC, CEM + LK, CEM with `blowing = LK` (336 states); the same grid with `boiling = ZGR` (29 of its states in the boiling regime) and one with `heat = NoHeat` in the regime (337) | every F finite; the surface solve brackets its root (no failed solve) | — | 673 states: 0 non-finite, 0 failed, 0 widened brackets | — (the never-NaN probe) |
| (f) | (d)'s state | informational: χ(T̄) 2.08391747, χ(T_s) 2.07981321 (−0.20 %), T_s − T̄ 3.68231670 K | — | — | — |
| (g) | (d)'s state through `rhsEvapBreakup`, Pilch–Erdman active (σ 5e-5 N/m, C_d 1, B 0.01: F(9)/Z(9) = 496.000253 1/s), μ_l 1e-3 as the breakup viscosity, `mu-liq` 0.1 (a decoy) | F(7) equals (d)'s model-2 value; `F(8) + Z(8) F(9)/Z(9)` equals (d)'s F(8) | 1e-14 rel.; 1e-13 rel. | 0; 4.2e-15 | the P2T branch dropped from model 4: 3.95e-2; `mu-liq` taken as μ_l: 4.48e-2; the carry paired with the total mass rate: 6.917e-2 |
| (h) | (d)'s state with TAB (an event model: model 2), μ_l 1e-3 as the breakup viscosity, `mu-liq` 0.1 | `F(1:8)` equal (d)'s | 1e-14 rel. | 0 | `mu-liq` taken whatever the breakup model: 1.918e-1 |
| (i) | `boiling = ZGR` at A = 0 (T_g = T_B = 600 K, `invTboil` stepped to it), Re 0, d 50 µm, k_l 0.6, ρ 958, c 4186; T̄ = T_B + 30, + 10, + 1, + 0.2 K | `TsOut` = T_B; F(8) = −10π k_l d ΔT_f/L_v (the flash fed by the parabolic conduction): −1.2510767205e-8, −4.1702557349e-9, −4.1702557349e-10, −2.0851278674e-10 kg/s (the 0.5 K floor at + 0.2 K); `(F(8)/F(7))/m = c/(L_v − c ΔT)` with the true ΔT | 1e-13 rel. | 4.4e-16; 2.2e-16 | the boiling heat-transfer coefficient's feed in place of the conduction: 9.947e-1 (F(8) at 5.3e-3 of its value at + 0.2 K); the vapour at the mean temperature: 5.557e-2; the floor dropped: 0.6 |
| (j) | (a)'s state with ε 0.8, through `rhsEvaporation` and through `rhsEvapBreakup` (Reitz–Diwakar at slip 0: F(9) = 0) | `TsOut` = the root of `10π k_l d (T_s − T̄) = π d k_g 2 (T_g − T_s) + ε σ π d² (T_g⁴ − T_s⁴)`, which the test solves by Newton (350.198043716 K); `F(7) = [π d k_g 2 (T_g − T_s) + ε σ π d² (T_g⁴ − T_s⁴)]/(m c)` | 1e-13 rel. | 2.2e-16, 2.2e-16 (both) | the radiative heat left out of the surface balance: T_s 5.655e-4, F(7) 7.905e-4; the radiation line of the right-hand side at T̄: F(7) 2.868e-4; that line deleted: F(7) 4.734e-3 |

Independent twins (python, the same states): T_s 323.682316697 K, F(7) 2.6403581310e4 K/s, h state 1.1036696988e8
J/kg/s, F(9)/Z(9) 496.000253 1/s (the surface balance by `brentq`, the A–S formulas typed); leg (j)'s root
350.19804371602882 K (mpmath, 40 digits).

The surface solve brackets the root exactly: `G(T̄) = −Q_net(T̄)/K` and `G(T_1) = (Q_net(T̄) − Q_net(T_1))/K` at
`T_1 = T̄ + Q_net(T̄)/K` (clipped to the property table's range and the gas temperature) have opposite signs whenever
`Q_net` does not increase with `T_s`; otherwise the bracket is widened to the gas temperature (heating) or the table's
lower end (cooling), and Brent's method runs to `4ε|T_s|`. A solve with no bracket returns a non-finite surface
temperature and the right-hand side a penalty slope, so the integrator rejects the step. Moving `T_1` onto `T̄` leaves
leg (a′) green: the widened bracket recovers the root.
